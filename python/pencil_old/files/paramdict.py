# Reads in F90 namelist as dictionary object. Optionally called by pc.read_param.
import re
r = re.compile(r'(?:[^,(]|\([^)]*\))+')

def param_formatter(part):
 part = re.sub(" ","",part)
 if part=="T": return True
 if part=="F": return False
 try:
  if "." in part: return float(part)
  else: return int(part)
 except: return re.sub("'","",part)

def tuplecatch(st):
 if "(" in st:
  st = st.replace("(","").replace(")","").split(",")
  for j in range(len(st)):
   st[j] = param_formatter(st[j])
  return tuple(st)
 else: return param_formatter(st)

def ReadNML(filename,nest=False):
 Params = dict()
 for rawline in open(filename):
  line = rawline.rstrip('\n')
  if line[0]=="&":
   supername=line[1:].lower()
   if nest: Params[supername] = dict()
  else:
   line = re.sub("^ ","",line)
   if line!="/":
    split = re.split("=",line)
    name = re.sub(" ","",split[0].lower())
    value = split[1]
    parts = r.findall(value)
    value = []
    for i in range(len(parts)):
     if "*" in parts[i]:
      s = parts[i].split("*")
      for i in range(int(s[0])):
       value += [tuplecatch(s[1])]
     else: value += [tuplecatch(parts[i])]
    if len(value)==1: value = value[0]
    if nest: Params[supername][name] = value
    else: Params[name] = value
 return Params
