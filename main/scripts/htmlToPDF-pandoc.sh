#!/bin/bash
filename=`basename $1 .html`
pandoc -t latex $filename.html \
	-V geometry:"top=2cm, bottom=1.5cm, left=1cm, right=1cm"\
    -V 'mainfont:DejaVuSerif.ttf' \
    -V 'sansfont:DejaVuSans.ttf' \
    -V 'monofont:DejaVuSansMono.ttf' \
    -V 'mathfont:texgyredejavu-math.otf' \
    -o $filename.pdf 
