# dockstore_sequenza

The workflow is made to run in Docker and uploaded to [Dockstore](https://docs.dockstore.org/en/develop/getting-started/getting-started.html).
You can find OICR's Dockstore page [here](https://dockstore.org/organizations/OICR).
The Docker container is based on [Modulator](https://gitlab.oicr.on.ca/ResearchIT/modulator), which builds environment modules to set up the docker runtime environment.

### Set Up and Run
Currently, this WDL must be run with Cromwell. 
It uses Cromwell configuration files to mount a directory to the docker container.
The directory contains data modules built with Modulator, which the WDL tasks need to access.
In addition, you must obtain run files locally and build data modules to a local directory.

#### 1. Build Data Modules
- Create a local `data_modules/` directory to store the data modules
    - make sure you have enough disk space
    - each data module could be 5-30 GB in size
- In future iterations of the workflow, this process will be simplified
- Enter the container:
```
# Mount this repository as /pipeline/; mount the data module destination directory as /data_modules/
docker run -it --rm -v [PWD]:/pipeline -v [data_modules]:/data_modules [CONTAINER ID (find in options.json)]

# Copy prerequisite code module YAMLs into the Modulator code directory (code/gsi/)
cp /pipeline/recipes/sequenza_data_modules_prep.yaml code/gsi/data_modules_recipe_prep.yaml
 
# Build the prerequisite code modules
./build-local-code code/gsi/data_modules_recipe_prep.yaml --output /data_modules --initsh /usr/share/modules/init/sh
 
# Copy data module YAMLs into the Modulator data directory (data/gsi/)
cp /pipeline/recipes/sequenza_data_modules.yaml data/gsi/data_modules_recipe.yaml
 
# Build the data modules
./build-local-data data/gsi/data_modules_recipe.yaml --output /data_modules --initsh /usr/share/modules/init/sh
 
# Change resulting file permissions
find /data_modules/ -type d -exec chmod 777 {} \; && \
find /data_modules/ -type f -exec chmod 777 {} \;
 
# /data_modules/ should now contain gsi/modulator/modulefiles/Ubuntu18.04/ and gsi/modulator/modulefiles/data/
```
For run directories that are not part of modules, copy them from UGE's archive at `/.mounts/labs/gsi/src/`

#### 2. Obtain Files Locally
In the test json, change file paths like so:
- File type files should be copied to local
    - E.g. use scp to copy from UGE
    - In the json, change the file path from UGE to local path
- String type files should be copied or moved to the mounted directory, if it's not already part of a module
    - In the json, change the file path to how the file would be accessed from inside the docker container
- $MODULE_ROOT paths can stay the same
```
# File type files
# File is copied to local machine
UGE:       "/.mounts/labs/gsi/testdata/wgsPipeline/input_data/wgsPipeline_test_pcsi/hg19_random.genome.sizes.bed"
Dockstore: "/home/ubuntu/data/sample_data/callability/hg19_random.genome.sizes.bed"
 
# String type files
# /data_modules/ is a directory mounted to the docker container
UGE:       "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.fa"
Dockstore: "/data_modules/gsi/modulator/sw/data/hg19-p13/hg19_random.fa"
 
# Root type paths
# The value of $MODULE_ROOT changes, but the path stays the same
UGE:       "$HG19_BWA_INDEX_ROOT/hg19_random.fa"
Dockstore: "$HG19_BWA_INDEX_ROOT/hg19_random.fa"
```

#### 3. Run with Cromwell
Submit the preprocessed subworkflow and modified json to Cromwell, with configs and options attached
```
# Validate the wrapper workflow and json
java -jar $womtool validate [WDL] --inputs [TEST JSON]
 
# For example:
java -jar $womtool validate wgsPipeline.wdl --inputs tests/wgsPipeline_test_cre_uge.json

# Submit to Cromwell
java -Dconfig.file=[CONFIG] -jar $cromwell run [WRAPPER WDL] --inputs [JSON] --options [OPTIONS]
 
# For example:
java -Dconfig.file=local.config -jar $cromwell run wgsPipeline.wdl --inputs tests/wgsPipeline_test_cre.json --options options.json
```