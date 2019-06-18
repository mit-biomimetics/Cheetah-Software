#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
from pyparsing import *
import hunspell
hobj = hunspell.HunSpell('/usr/share/hunspell/en_US.dic', '/usr/share/hunspell/en_US.aff')
import subprocess
from multiprocessing import Process, Queue, Lock, Pool, Manager
from itertools import chain

lt_exclusions = [
  'WHITESPACE_RULE',
  'EN_A_VS_AN',
  'COMMA_PARENTHESIS_WHITESPACE',
  'GENERAL_XX',
  'SENTENCE_WHITESPACE',
  'EN_QUOTES',
  'DT_DT',
  'A_PLURAL',
  'DOUBLE_PUNCTUATION',
  'ALL_OF_THE',
  'ENGLISH_WORD_REPEAT_BEGINNING_RULE',
  'UPPERCASE_SENTENCE_START',
  'ENGLISH_WORD_REPEAT_RULE',
  'IN_A_X_MANNER',
  'HE_VERB_AGR',
  'PHRASE_REPETITION',
  'KIND_OF_A',
  'EN_UNPAIRED_BRACKETS',
  'MAY_BE'
]

import sys
import re
import os

def checkfile(f,lock):
  localsuccess = True
  text = ""
  spells = ""
  interface = None
  m = re.search("interfaces/(\w+)/",f)
  if m:
    interface = re.compile(r"\b"+m.group(1)+r"\b",re.I)
  prevs = ("",0,0)
  fc = file(f,'r').read()
  if f.endswith(".tex"):
    chunks = iter([(fc,0,0)])
  else:
    chunks = chain(cppStyleComment.scanString(fc),dblQuotedString.scanString(fc))
  for (m,s,e) in chunks:

    
    lineno = len(fc[:s].split("\n"))
    t = m if isinstance(m,str)  else m.asList()[0]
    #if re.match("^\s+$",fc[prevs[2]:s]) and not t.startswith('///') and t.startswith('//') and prevs[0].startswith('///') and "@{" not in prevs[0] and "@}" not in prevs[0]:
    #  print "%s:%d %s" % (f,lineno,m)
    #prevs = (t,s,e)
    if t[1:-1].endswith(".hpp") or t[1:-1].endswith(".h"):
      continue
      
    
    if t.startswith("/*") and not t.startswith("/**"):
      continue
    #if "//" in t[3:-3]:
    #  continue
    
    t = re.sub(r'\bNOLINT\(.*?\)','',t)

    t = re.sub(r'\s[^\s]+(-by-)[^\s\)]+(\s)','\\2',t)
    
    #tex
    t = re.sub('\\\\begin\{(equation|verbatim|pytex)\}.*?\\\\end\{\\1\}','',t,flags=re.DOTALL)
    t = re.sub('\\\\(newcounter|setcounter|addtocounter|texttt|arabic|alph|ref|label|textit|emph)\{.*?\}','foo',t,flags=re.DOTALL)
    t = re.sub('\\\\verb\|.*?\|','',t,flags=re.DOTALL)
    t = re.sub('[\d\.]+\w+','',t)
    t = re.sub('\\\\begin\{(.*?)\}(.*?)\\\\end\{\\1\}','\\2',t,flags=re.DOTALL)
    t = re.sub('\\\\(usepackage|includegraphics).*?\{.*?\}','',t,flags=re.DOTALL)
    t = re.sub('\$.*?\$','',t)
    
    t = re.sub('\[\w+\]','',t)
    t = re.sub('\(\w+\s+x\s+\w+\)','',t)
    t = re.sub('[\\\\@](copydoc|a|e|p|param|defgroup)\s+\w+','',t)
    t = re.sub(r'\b([A-Z]-)[A-Z_]{2,}\w+','',t)
    t = re.sub(r'(?<!\w)-+(?!\w)','',t)
    t = re.sub('\\\\verbatim(.*?)\\\\endverbatim','',t,flags=re.DOTALL)
    t = re.sub('\\\\f\$(.*?)\\\\f\$','',t,flags=re.DOTALL)
    t = re.sub('\\\\f\[(.*?)\\\\f\]\$','',t,flags=re.DOTALL)
    

    
    t = re.sub(r'\b\w+[A-Z]\w+\b','',t) # camelcase
    t = re.sub(r'\b\w+_\w*\b','',t) # camelcase
    t = re.sub('[\\\\@][a-zA-Z]+\{.*?\}','',t)
    t = re.sub('[\\\\@][a-zA-Z{}]+','',t)
    
    t = re.sub('<tt>.*?</tt>','',t)
    t = re.sub('</?\w+>','',t)
    t = re.sub('#\w+','',t)
    t = re.sub('\*','',t)

    t = re.sub('C\d+','',t) # Visual studio warnings
    t = re.sub("'[^\s]*?'",'',t)
    t = re.sub(r'\b[ntdpcx][a-zA-Z0-9]\b',lambda e: e.group(0) if hobj.spell(e.group()) else '',t) # should be escaped with \e
    t = re.sub("-\d+",'',t)
    t = re.sub("\b\w+\d+\b",'',t)
    t = re.sub("- ",' ',t)
    if interface is not None:
      t = interface.sub('',t)
      
    
    text+= t+"\n"
    #for w in t.split(" "):
    #  if not hobj.spell(w):
    #    success = False
    #    spells+="Spelling error: %s (%s)\n" % (w,str(hobj.suggest(w)))
    prevs = (t,s,e)
    
  p = subprocess.Popen(['java','-jar','/home/jg/programs/LanguageTool-2.5/languagetool-commandline.jar','-l','en','-d',",".join(lt_exclusions),'-'],stdin = subprocess.PIPE, stdout = subprocess.PIPE)
  out, err = p.communicate(text)
  if len(out.split("\n")[2:-2]) > 0:
    localsuccess = False
  ph = subprocess.Popen(['hunspell','-l','-p','casadi.dic'],stdin = subprocess.PIPE, stdout = subprocess.PIPE)
  outh, errh = ph.communicate(text)
  
  if len(outh) > 0:
    localsuccess = False
  if not localsuccess:
    lock.acquire()
    print "In file %s" % f
    print "\n".join(out.split("\n")[2:-2])
    print outh
    lock.release()
    return False
  else:
    return True
  

if __name__ == "__main__":
  pool = Pool(processes=16)

  dir = sys.argv[1]
  files = sys.argv[2:]

  if 'core' in dir:
    sys.exit(0)

  success = True
  
  results = []

  queue = Queue()

  manager = Manager()
  lock = manager.Lock()

  for f in files:
    localsuccess = True
    ff = os.path.join(dir,f)
    if f.endswith("cpp"):
      continue
    results.append(pool.apply_async(checkfile, args=(ff,lock)))
    
  pool.close()
  pool.join()
  
  success = True
  
  for r in results:
    if not r.get():
      success = False

  if success:
    sys.exit(0)
  else:
    sys.exit(1)
