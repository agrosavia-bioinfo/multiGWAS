#!/bin/bash
filename=`basename $1 .html`
pandoc -s $filename.html -o $filename.docx
