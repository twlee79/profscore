from docutils.core import publish_string
infile = open("neuroserpin/neuroserpin_example.rst")
input = infile.read()
input = input.replace('example/neuroserpin/','')
outfile1 = open("neuroserpin/neuroserpin_example.txt","w")
outfile1.write(input)
outfile2 = open("neuroserpin/neuroserpin_example.html","w")
output = publish_string(input,writer_name='html')
#, settings_overrides={'stylesheet':None,'stylesheet_path':'html4css1.css'}
outfile2.write(output)


infile = open("antithrombin/antithrombin_example.rst")
input = infile.read()
input = input.replace('example/antithrombin/','')
outfile1 = open("antithrombin/antithrombin_example.txt","w")
outfile1.write(input)
outfile2 = open("antithrombin/antithrombin_example.html","w")
output = publish_string(input,writer_name='html')
#, settings_overrides={'stylesheet':None,'stylesheet_path':'html4css1.css'}
outfile2.write(output)