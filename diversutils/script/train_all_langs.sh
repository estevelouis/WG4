langs=(eo en)
for i in "${!langs[@]}"
do
	lang="${langs[$i]}"
	echo "trying lang=$lang"
	bash script/train_vectors.sh $lang
done
