# Copyright (c) 2022 Michiel Bongaerts.
# All rights reserved.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.


"""
@author: Michiel Bongaerts
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures



class WeightedMeanSigma:
    def __init__(self,a,b):
        self.a = a
        self.b = b
    
    def fit( self,age_train, sigma_train ):
        assert( len( age_train ) == len( sigma_train ))
        self.age_train = age_train
        self.sigma_train = sigma_train
        
    def predict(self, ages):
        means = []
        for age in ages:
            weigths = np.exp( - np.abs( self.age_train - age )/ (self.a + self.b * age ))
            weigths = np.round(weigths,10)
            weigths = weigths / np.sum(weigths)
    
            weighted_mean = np.sum( (weigths * self.sigma_train)  )
            means.append(weighted_mean)
            
        return means    
    

class LR_age_sex:
    ''' For details about this model see https://doi.org/10.3390/metabo11010008 '''
    def __init__(self, n_jobs=1, polynomials = 3,a= 35,b= 0.1, Z_outlier_threshold = 3 ):
        self.n_jobs = n_jobs
        self.LR = LinearRegression(n_jobs= self.n_jobs, fit_intercept=False)
        self.polynomials = polynomials
        self.a = a 
        self.b = b
        self.fitted = False
        self.Z_outlier_threshold = Z_outlier_threshold
        
    def fit(self, X_train, y_train, interaction = False):
        # First columnn of X_train needs to be age in days
        assert( X_train.columns.tolist()[0] == 'age_in_days' )
        
        self.fitted = True
        self.interaction = interaction
      
        self.y_train_median = y_train.median()
        self.y_train_mad = y_train.mad()
       
        # Determine outliers 
        y_train_Z = (y_train - self.y_train_median )/ ( self.y_train_mad * 1.48 )
        self.not_outlier_IDs = y_train_Z[ y_train_Z.abs() < self.Z_outlier_threshold ].index.tolist()

        # Add polynomials
        poly = PolynomialFeatures( self.polynomials ,include_bias=True)
        X_poly = pd.DataFrame( poly.fit_transform( X_train[['age_in_days']]),index= X_train.index )

        self.feature_names = X_poly.columns.tolist()


        # Remove outliers and use Z-tranformed y_train
        self.X_train_not_outliers = X_train.loc[self.not_outlier_IDs]
        self.X_train_poly = X_poly.loc[self.not_outlier_IDs]
        self.y_train = y_train.loc[self.not_outlier_IDs]
        self.ages = X_train['age_in_days'].loc[self.not_outlier_IDs].values


        self.LR.fit(self.X_train_poly.values, self.y_train.values)
        self.y_fit = self.LR.predict(self.X_train_poly.values)
        self.sigma_train = (self.y_train - self.y_fit)**2
        
        self.WMS = WeightedMeanSigma(self.a, self.b)
        self.WMS.fit(self.ages,self.sigma_train)

        self.sigma_WMS = self.WMS.predict( self.ages ) 

        self.Sigma = np.diag( self.sigma_WMS )
        self.XTX_inv = np.linalg.inv(np.dot(np.transpose(self.X_train_poly.values) , self.X_train_poly.values))
        self.cov_beta = np.dot( np.dot( np.dot( np.dot(self.XTX_inv, self.X_train_poly.values.T), self.Sigma ), self.X_train_poly.values ), self.XTX_inv )


        
    def predict(self, X_test):
        # First columnn of X_train needs to be age in days
        assert( X_test.columns.tolist()[0] == 'age_in_days' )
    
        # prevent extrapolating on age which much higher than ages included in fit
        max_age = self.X_train_not_outliers['age_in_days'].max()
        X_test.loc[ X_test['age_in_days'] > max_age, 'age_in_days' ]  = max_age
    
        poly = PolynomialFeatures( self.polynomials ,include_bias = True)
        X_poly = pd.DataFrame( poly.fit_transform( X_test[['age_in_days']]),index= X_test.index )


        self.X_test = X_poly.values
        self.ages_test = X_test['age_in_days']
        
        self.y_pred = self.LR.predict(self.X_test) 


        self.var_1 = np.array( [np.dot(np.dot(g, self.cov_beta), g) for g in self.X_test] )
        self.var_2 = self.WMS.predict(self.ages_test) 
        self.std_pred = np.sqrt(  self.var_1 + self.var_2  ) 

        results = pd.DataFrame(self.y_pred , columns=['y_pred'])
        results.loc[:,"std_pred"] = self.std_pred


        return results
        

class OligoPipeline:
    
    # distionary with all compounds per disease
    disease_with_compounds = {  'GM1' : ['945_GM1', '909_GM1','1472_GM1','1436_GM1','1071_GM1','1274_GM1'],
                                'GM2' : [ '783_GM2','747_GM2','1351_GM2','1315_GM2','1148_GM2','1112_GM2' ],
                                'α-mannosidosis': ['1030_α-mann','904_α-mann','868_α-mann','580_α-mann','544_α-mann'],
                                'α-NAGA': ['321_α-NAG','307_α-NAG' ], 
                                'Sialic acid storage disease' : ['N-acetylneuraminic acid'],
                                'GSD II' : ['1151_Hex7','665_Glc4','989_Hex6'],
                                'AGU' : ['787_aspartylglucos','496_aspartylglucos','334_aspartylglucos'],
                                'Fucosidosis' : ['480_fuco','366_fuco','1055_fuco'],
                                'β-mannosidosis' : ['418_β-mann','382_β-mann'],
                                'Galactosialidosis' : ['1200_(galacto)sia','1008_(galacto)sia'],
                           
                                }
    
    compound_with_disease = { c:d for d,cs in disease_with_compounds.items() for c in cs  }
        
    
    # Fix the order of compounds
    compound_disease_sorted = [ 'N-acetylneuraminic acid','1351_GM2','1315_GM2','783_GM2','747_GM2',
                                '1148_GM2','1112_GM2','1472_GM1','1436_GM1','1274_GM1','1071_GM1','945_GM1','909_GM1','1055_fuco','480_fuco',
                                '366_fuco','418_β-mann','382_β-mann','787_aspartylglucos','496_aspartylglucos',
                                '334_aspartylglucos','321_α-NAG','307_α-NAG','1030_α-mann','904_α-mann','868_α-mann', 
                                '580_α-mann','544_α-mann','1151_Hex7','989_Hex6','665_Glc4','1200_(galacto)sia','1008_(galacto)sia']

    internal_standards = ['1235_Hexaacetyl-chitohexaose', '1271_hexaacetyl-chitohexaose', '644_acarbose', '680_acarbose','13C3 N-acetylneuraminic acid']
    
    def __init__(self, data_batch, all_ref_data):
        
        self.check_input_data(data_batch)
        self.check_input_data(all_ref_data, ref=True)
        
        self.data_batch = data_batch
        self.all_ref_data = all_ref_data
        self.reference_values_determined = False
        
    def check_input_data(self, df, ref=False):
        
        columns = set(['compound', 'sample_amt', 'creatinine_mmol_liter', 
                       'age_in_years','sample_ID', 'sample_amt_normalized', 
                       ])
        
        if( ref == True):
            columns = columns.union(['batch_name', 'group'])
        
        overlap = columns.difference(df.columns)
        
        if( len(overlap) > 0):
            print(overlap)
            raise AssertionError('Input dataframe does not contain the right column(s)(names)')
            


    
    def determine_reference_values(self, plot = True):
        
        all_ref_data = self.all_ref_data
        
        # Make plot where we see the regression model fitted for establishing reference values
        if( plot == True):
            fig = plt.figure(figsize=(40,40))
            fig.subplots_adjust(hspace=0.4, wspace=0.5)


        self.reference_values = {}
        for i,(compound, data_gb) in enumerate( all_ref_data.loc[ (all_ref_data['group'] == 'control') & (all_ref_data['compound'].isin(OligoPipeline.compound_disease_sorted) )].groupby('compound') ):
            
            data_gb = data_gb.assign(age_in_days = data_gb['age_in_years'] *365)
            data_gb = data_gb.sort_values(by='age_in_years')
            data_gb = data_gb.loc[ data_gb['sample_amt_normalized'].notnull()]

            # Removal of outliers
            Z = data_gb['sample_amt_normalized']
            Z = ( Z- Z.median()) / (Z.mad() )
            data_gb = data_gb.loc[ Z.abs() < 5 ]
    
    
            if( plot == True):
                ax = fig.add_subplot(8,8,i+1)

            # If we have enough non-zero numbers we continue with fitting the regression model
            if( sum(data_gb['sample_amt_normalized'] > 0)  >= 100 ):


                X = data_gb.set_index('sample_ID')[['age_in_days']]
                y = data_gb.set_index('sample_ID')['sample_amt_normalized']

                # Fit regression model
                lr = LR_age_sex(a=100, b=1, Z_outlier_threshold=3)
                lr.fit(X[['age_in_days']],y)

                # Make new dataframe for different ages to predict the mean and std from regression model
                X_pred = np.linspace( X['age_in_days'].min(), X['age_in_days'].max(),100)[:,np.newaxis]
                X_pred = pd.DataFrame(X_pred.tolist(), columns = ['age_in_days'])
                y_pred = lr.predict(X_pred)

                if( plot == True):
                    # Plot predictions
                    ax.plot( X_pred['age_in_days']/365, y_pred['y_pred'] ,color='red',alpha=1,zorder=1000,linewidth=4 )
                    ax.plot( X_pred['age_in_days']/365, y_pred['y_pred'] + y_pred['std_pred'] ,color='red',alpha=1,zorder=1000,linewidth=2 )
                    ax.plot( X_pred['age_in_days']/365, y_pred['y_pred'] - y_pred['std_pred'] ,color='red',alpha=1,zorder=1000,linewidth=2 )

                    ax.fill_between(X_pred['age_in_days']/365, 
                                    y_pred['y_pred'] - y_pred['std_pred'],
                                    y_pred['y_pred'] +y_pred['std_pred'], alpha=0.1,color='darkred')


                self.reference_values[compound] = lr

            # Not enough data for fitting for regression so just take 'naive' std and mean
            else:
                mean =  data_gb['sample_amt_normalized'].mean()
                std = data_gb['sample_amt_normalized'].std() + 0.1
                self.reference_values[compound] = {'mean': mean,
                                       'std':  std 
                                      }
                
                if( plot == True):
                    ax.hlines(mean,0,70, color='red',linewidth=3)
                    ax.hlines(mean-std,0,70, color='red',linewidth=3)
                    ax.hlines(mean+std,0,70, color='red',linewidth=3)

            if( plot == True):
                ax.scatter(  data_gb['age_in_days']/365, data_gb['sample_amt_normalized'] )
                ax.set_title(compound+" N={}".format(data_gb.shape[0]))
                ax.set_ylabel('Sample amount \n normalized \n by creatinine',fontsize=10 )
                ax.set_xlabel('Age (years)',fontsize=10 )
                ax.set_xlim([-1,70])
                
        self.reference_values_determined = True

            
    def determine_Z_scores(self, df ):
        assert( self.reference_values_determined == True)
        
        def Z_scores_ref(row):
            compound = row['compound']

            if( compound in self.reference_values):

                # Check if regression was made on this compounds
                if( isinstance( self.reference_values[compound], LR_age_sex) ):
                    ID = row['sample_ID']
                    age_in_days = row['age_in_years'] * 365

                    lr = self.reference_values[compound]

                    X = pd.DataFrame([[age_in_days]],columns= ['age_in_days'])
                    y_pred = lr.predict( X )

                    mean = y_pred['y_pred'].tolist()[0]
                    if( mean < 0): mean = 0
                    std = y_pred['std_pred'].tolist()[0]

                    Z = (row['sample_amt_normalized'] - mean) / std
                    LR_or_fixed = 'LR'


                # Else take other reference (fixed) reference values
                else:
                    mean = self.reference_values[compound]['mean']
                    if( mean < 0): mean = 0
                    std = self.reference_values[compound]['std']
                    if( std <= 0): std = 0.1

                    Z = (row['sample_amt_normalized'] - mean) / std
                    LR_or_fixed = 'fixed'

            else:
                Z =  np.nan
                mean = 0
                std = 0
                LR_or_fixed = 'No ref'

            return pd.Series({'reference_values': (np.round(mean,2),np.round(std,2), LR_or_fixed ), 'Z_scores': Z})

    
        df = df.assign(**df.apply(lambda row: Z_scores_ref( row[['compound','sample_amt_normalized','sample_ID','age_in_years']] ), axis=1))
        return df
    
    
    

    def plot_Z_scores_patient(self,df, sample_ID):
        assert( self.reference_values_determined == True)
        assert( 'Z_scores' in df.columns )
        
        diseases = list(OligoPipeline.disease_with_compounds.keys())
        colors = sns.color_palette('rainbow',len(diseases)*10)
        colors = np.array( colors )
        colors = colors[np.arange(0,len(colors),10)]

        # make a dictionary with a color for every disease
        diseases = np.unique( list( OligoPipeline.compound_with_disease.values() ))

        disease_with_colors = { u:colors[-i-1] for i,u in enumerate(diseases) }
        compound_with_color = { c:disease_with_colors[d] for c,d in OligoPipeline.compound_with_disease.items() }


        # This could be used to scale the Z-scores but when you don't want any scaling then just return x itself
        def Z_score_transform(x):
            return np.sign(x) * np.log(abs(x) + 0.5)


        data_gb = df.loc[ df['sample_ID'] == sample_ID]

        sns.set(font_scale=1,style='white', rc = {'font.family' : 'serif'})
        fig = plt.figure(figsize=(20,15))
        ax = fig.add_subplot(1,2,1)

        # Apply transformation on Z-scores
        data_gb = data_gb.assign(Z_scores_scaled = data_gb['Z_scores'].apply( lambda x: Z_score_transform(x)  ) )


        # Check if all compounds are present. If not, add compound with NaNs
        missing_compounds = list( set(OligoPipeline.compound_disease_sorted).difference( set(data_gb['compound'] ) ) )
        if( len(missing_compounds) > 0):
            print("There a missing compounds for {}, which are: {}".format(sample_ID, missing_compounds) )

            row = { col : np.nan for col in data_gb.columns }
            for c in missing_compounds:
                row['compound'] = c
                # add new row in dataframe to have at least the compound with NaNs in there.
                data_gb = data_gb.append( pd.Series(row), ignore_index = True )


        # plot values current sample
        data_gb = data_gb.set_index('compound').loc[OligoPipeline.compound_disease_sorted]
        compound_colors = [disease_with_colors[ OligoPipeline.compound_with_disease[el]  ] for el in OligoPipeline.compound_disease_sorted ]

        # Put a minus in front of position are correctly displayed
        positions = [ pos for pos in range(len(data_gb['Z_scores_scaled'])) ]

        # makes horizonal lines from zero axis to Z-score
        for x,y,color,compound_disease in zip(data_gb['Z_scores_scaled'], positions, compound_colors, OligoPipeline.compound_disease_sorted ):
            ax.plot( [0,x], [y,y], c=color, zorder= -1000, linewidth=4)



        # plot the Z-scores as squares 
        ax.scatter( data_gb['Z_scores_scaled'], positions , c = compound_colors, 
                   edgecolor='k',linewidth=2,s=120,marker='s',zorder= 1000)

        # Plot names of compounds
        ax.set_yticks( positions )
        ax.set_yticklabels( data_gb.index.tolist() )

        # obtain current x,y limits for plot
        ylims = ax.get_ylim()
        xlims = ax.get_xlim()

        # plot separation lines for diseases
        disease_seperation_line_plotted = []
        for x,y,color,compound_disease in zip(data_gb['Z_scores_scaled'], positions, compound_colors, OligoPipeline.compound_disease_sorted ):

            # Plot a line to separate the disease markers
            disease = OligoPipeline.compound_with_disease[compound_disease]

            # Check if line is already plotted if not plot 
            if( disease not in disease_seperation_line_plotted):
                ax.hlines( y - 0.5, -10**10, 10**10 ,color='lightgrey')
                disease_seperation_line_plotted.append(disease)

            ax.set_title(sample_ID, fontsize=20, y=1.01)

        ax.set_xlabel('Z-score',fontsize=25)


        # create some lines
        zs = [1,10,100,1000,10000,100000,1000000]
        ax.vlines( 0, -10**10, 10**10, color='lightgrey',zorder=-1000 )
        for z in zs:
            ax.vlines(   Z_score_transform(z),-10**10,10**10, color='lightgrey',zorder=-1000 )
            ax.vlines(  -Z_score_transform(z),-10**10,10**10, color='lightgrey',zorder=-1000 )

        x_ticks = [ Z_score_transform(-z) for z in zs ]
        x_ticks.append(0)
        x_ticks.extend(  [ Z_score_transform(z) for z in zs ] )
        ax.set_xticks( x_ticks )

        x_ticks_labels = [ str(-z) for z in zs ]
        x_ticks_labels.append('0')
        x_ticks_labels.extend( [ str(z) for z in zs ] )
        ax.set_xticklabels( x_ticks_labels )


        ax.set_xlim([ Z_score_transform(-10), max(xlims) *1.1 ])
        ax.set_ylim(ylims)



        # Plot values of patient with reference values
        ax = fig.add_subplot(1,2, 2)

        ax.set_yticks( positions )
        ax.set_yticklabels( ['' for i in range(len(positions))]) 
        y_ticks = ax.get_yticks()
        ax.set_ylim( ylims )

        def clean(x):
            chars = ['(',')',"'",' ']
            for char in chars:
                x = x.replace(char,'') 
            return x

        a, b,c,d = ['Patient', 'Mean',   'Std',  'Model']

        for ii,el in enumerate([a,b,c,d]):
            ax.text(ii*0.15, max(y_ticks)+1,'{}'.format(el),weight ='bold')

        last_color = None
        for y,ref,A_pat,color in zip(positions,
                                     data_gb['reference_values'], 
                                     data_gb['sample_amt_normalized'],
                                     compound_colors):        

            if( str(ref) == 'nan'):
                ref = 'None, None, None'

            a = round(A_pat,2)
            b,c,d = [ clean(el) for el in str(ref).split(',')]
            d = d[0].capitalize() + d[1:]
            if(d == 'LR'):
                d = 'Regression'

            for ii,el in enumerate([a,b,c,d]):
                ax.text(ii*0.15,
                        y,
                        '{}'.format(el),
                        backgroundcolor='w',weight ='normal',zorder=ii,)

            last_color = color


        x_l,x_u = ax.get_xlim()
        x_u = ii*0.15 + 0.15
        # plot separation lines for diseases
        disease_seperation_line_plotted = []
        for x,y,color,compound_disease in zip(data_gb['Z_scores_scaled'], positions, compound_colors, OligoPipeline.compound_disease_sorted ):

            # Plot a line to separate the disease markers
            disease = OligoPipeline.compound_with_disease[compound_disease]

            # Check if line is already plotted if not plot 
            if( disease not in disease_seperation_line_plotted):
                ax.hlines( y - 0.5, x_l,x_u- 0.05 ,color='lightgrey')
                disease_seperation_line_plotted.append(disease)

        ax.set_xlim(x_l,x_u)


        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)
        ax.set_frame_on(False)

        return fig
    