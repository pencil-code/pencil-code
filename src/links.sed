/^!\*\* AUTOMATIC REFERENCE-LINK/, /^\*\*\*\*/ {
#/^! *.*, %/ { 
s/^! *\(.*\), *%.*$/https:\/\/ui.adsabs.harvard.edu\/abs\/\1/ 
t end
d
b cont
: end
p
d
: cont
}
d
