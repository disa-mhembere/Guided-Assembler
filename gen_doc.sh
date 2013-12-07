#!/bin/bash

# gen_doc.sh
# Email: disa@jhu.edu, sleee320@jhu.edu
# Copyright (c) 2013. All rights reserved.

# Generate awesome docs & clean files we don't want

epydoc -v -o doc/ --pdf -n Guided\ Assember\ API src/*.py tests/*.py --debug
cd doc/ && find . | grep -v "api.pdf" | xargs rm
