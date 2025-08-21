# This project is used to gain familiarity with Docker and Nextflow with Python


# Installation (working in VS Code and WSL within Windows 11)

1. Project directory was made within WSL partition: 
```
/home/oneohtrix/projects/demo_workflow
```

2. I had to install miniconda within WSL
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

3. Specified a YAML file to set up the environment (in WSL) and initiated
```
/home/oneohtrix/projects/demo_workflow/environment.yml
conda env create -f environment.yml
```

4. Initiated git and .gitignore
```
git init
touch .gitignore
```

5. Needed to install WSL extension in VS Code, otherwise IDE runs in Windows 11
and doesn't recognize the installed environment+python. VS Code then can connect
to the WSL as a server by opening terminal and typing:
```
code .
```

6. Acquired a Google AI API key. Note that this requires setting up a Google 
account distinct from Google One/Pro
```
https://aistudio.google.com/apikey
AIzaSyCz9h2gUpSSqy-_JFhHOurD6WLeO-wF9Hk
```
# aggregate_sc_coexpression_wdl
# aggregate_sc_coexpression_wdl
