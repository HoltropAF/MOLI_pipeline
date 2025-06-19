# Quick Conda + VS Code Setup

_Assumes Python and Anaconda/Miniconda are installed and available in your system's PATH._

## Set Up Conda Environment

1. Place `conda_oligo_github.yml` (or relevant `_exe_vX` file) in your project folder  
2. In terminal:  
   `cd /path/to/folder`
3. Create environment:  
   `conda env create -f conda_oligo_github.yml`  
   (Optional custom name: `-n my_env`)
4. Activate:  
   `conda activate <env_name>`
5. In VS Code:  
   `Ctrl+Shift+P` → `Python: Select Interpreter` → choose your environment  
   (More info: [VS Code Python environments](https://code.visualstudio.com/docs/python/environments))
6. ✅ Done

_Note: If the environment doesn’t appear in VS Code, restart VS Code or install the Python extension._

For more on creating environments from YAML files:  
[Conda: Creating an environment from a .yml file](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)

---

## Pipreqs (for generating `requirements.txt`)

Used to log dependencies if changes are made to the environment.

1. Open Anaconda Prompt  
2. `cd` to the folder with your `.py` file  
3. `conda activate <env_name>`  
4. If needed: `pip install pipreqs`  
5. Run:  
   - For a single script: `pipreqs .`  
   - For a folder: `pipreqs . --force`  
   - Save elsewhere: `pipreqs /your/folder/path`

✅ Success message:  
`INFO: Successfully saved requirements file in ./requirements.txt`
