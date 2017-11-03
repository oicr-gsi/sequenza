#!/bin/bash
cd $1
ls | sed 's/.*\.//' | sort | uniq -c
