#!/usr/bin/python
import sys, os

if len(sys.argv)<2:
   sys.stderr.write('Usage: %s CLASS_NAME\n\n' % sys.argv[0])
   sys.exit(1)
if not 'MRB_TOP' in os.environ.keys():
   sys.stderr.write('$MRB_TOP not defined!\n\n')
   sys.exit(1)

name='CMAlgo' + sys.argv[1]
working_package='CMTool'
target_dir='%s/srcs/larreco/RecoAlg/CMTool' % os.environ['MRB_TOP']
source_dir='%s/mac/tmp' % target_dir

in_source='%s/cmalgo.cxx' % source_dir
in_header='%s/cmalgo.h' % source_dir

src_list = { in_source   : '%s/%s.cxx'         % (target_dir,name),
             in_header   : '%s/%s.h'         % (target_dir,name) }

for src in src_list.keys():
   if os.path.isfile(src_list[src]):
      sys.stderr.write('File already exists: %s\n\n' % src_list[src])
      sys.exit(1)
      
for src in src_list.keys():
   contents=open(src,'r').read()
   contents=contents.replace('CMALGO_CLASS_NAME',name.upper())
   contents=contents.replace('cmalgo_class_name',name.lower())
   contents=contents.replace('CMAlgo_Class_Name',name)
   contents=contents.replace('USER',os.environ['USER'])
   contents=contents.replace('Working_Package',working_package)
   fout=open(src_list[src],'w')
   fout.write(contents)
   fout.close()

print
print 'Generated the followings under %s.' % target_dir
for key in src_list.keys():
   print '    %s' % src_list[key]
print
sys.exit(0)
