#!/usr/bin/env julia
#
# Reads the Julia source files and makes an HTML manual out of
# the Markdown documentation.
# 
# Author: Daniel Carrera (danielc@astro.lu.se)
# 
# REQUIREMENTS:  python-markdown library
#                --> apt-get install python-markdown
# 
# Usage:     julia make_manual.jl

using PyCall
@pyimport markdown


text  = ""

md_files = [ "README.md" ]
jl_files = [ "src/timeseries.jl", "src/particles.jl", "src/averages.jl", "src/dimensions.jl" ]

for file = md_files
	info("Processing $file")
	
	src  = open(file, "r")
	text = text * readall(src)
	close(src)
end

for file = jl_files
	info("Processing $file")
	
	src  = open(file, "r")
	save = false
	for line = readlines(src)
		
		if beginswith(line,"#=doc")
			save = true
		elseif beginswith(line,"#=end")
			save = false
		elseif save
			if beginswith(line,"#")  start = 2  end
			if beginswith(line,"# ") start = 3  end
			text = text * line[start:end]
		end
		
	end
	close(src)
end

#
# Write the HTML documentation as one file.
#
header = "
<html>
<head>
  <link href='manual.css' rel='stylesheet'  type='text/css'>
</head>
<body>
";
footer = "
</body>
</html>
";

ext  = ["tables","fenced_code", "codehilite", "toc", "def_list", "extra"]
html = header * markdown.markdown(text, extensions=ext) * footer

#
# Alternate row styles.
#
count = 0
lines = split(html,"\n")
for i = 1:length(lines)
	if contains(lines[i],"<table>")
		count = 1
	end
	if contains(lines[i],"<tr>")
		style = count % 2 == 1 ? "odd" : "even"
		lines[i] = replace(lines[i],"<tr>","<tr class='$style'>")
		count += 1
	end
end

html = join(lines,"\n")
html = replace(html, r"<table>", "<table cellspacing='0' cellpadding='0'>")


doc = open("doc/manual.html", "w")
write(doc, html)
close(doc)
	
