prec.zip:
	zip -r prec.zip data/precalculated/
upload_prec:prec.zip
	scp ./prec.zip sunhf@compbio.tongji.edu.cn:~/public_html
download_prec:
	wget compbio.tongji.edu.cn/~sunhf/prec.zip
	unzip -o prec.zip
	
data.zip:
	zip -r data.zip data/AGP/ data/RNAseq/ data/immune_datasets/ data/test/
upload_data:data.zip
	scp ./data.zip sunhf@compbio.tongji.edu.cn:~/public_html
download_data:
	wget compbio.tongji.edu.cn/~sunhf/data.zip
	unzip -o data.zip

select:
	Rscript selection.R	

clean:
	rm *.zip

