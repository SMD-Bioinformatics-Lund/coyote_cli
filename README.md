## Import coyote data

### Functions
Load genetic data into coyote. This can be achieved in two ways.
* Using flags via the load command. This is similar to the legacy perl-script
* Using YAML-file via the yaml command.

This new CLI is written in python and combines all various data types that 4 different perl-scripts handled before. It also handles both DNA and RNA samples. Additionally improved case_id incrementation that was crudely handled by perl-scripts before and adds the option to update cases with new information or overwriting old information. This could be both meta-information or types of genetic variation.

Handles data types:
* CSV (gzipped)
* YAML
* JSON
* VCF

import_coyote_sample.py --help

import_coyote_sample.py load --help

import_coyote_sample.py yaml --help


### todo
* make sure update flag works for yaml import
* documentation
* better testing examples 
  * (SNV-indels that are not all filtered)
  * fewer of other data types
* create pip requirements.txt

	
