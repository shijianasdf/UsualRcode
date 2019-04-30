library(rmarkdown)
render(input='/pub6/temp/xxzhang/test/test.Rmd', output_format = 'html_document', output_file = 'testnew.html', output_dir = '/pub6/temp/xxzhang/test/')
render(input='/pub6/temp/xxzhang/test/test.Rmd', output_format = 'word_document', output_file = 'testnew.doc', output_dir = '/pub6/temp/xxzhang/test/')
render(input='/pub6/temp/xxzhang/test/test.Rmd', output_format = pdf_document(latex_engine = "xelatex", pandoc_args="--metadata=mainfont:SimSun", template="/pub5/xiaoyun/Bin/template.latex"), output_file = 'testnew.pdf', output_dir = '/pub6/temp/xxzhang/test/', encoding="UTF-8")
