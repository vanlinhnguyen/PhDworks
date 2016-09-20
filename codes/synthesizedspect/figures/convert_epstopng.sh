for f in `ls *.eps`; do
	convert -density 600 $f -flatten ${f%.*}.png;
done 


