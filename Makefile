.PHONY: pack
pack:
	zip -r data.zip data/precalculated/
raw:
	zip -r raw.zip data/AGP/ data/RNAseq/ data/immune_datasets/
upload:pack
	scp ./data.zip sunhf@compbio.tongji.edu.cn:~/public_html
clean:
	rm data.zip
