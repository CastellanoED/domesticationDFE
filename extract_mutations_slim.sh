FILE=$1
INIT=`cat $FILE | grep -n 'Mutations' | sed 's/:.*//g'`
END=`cat $FILE | grep -n 'Genomes' | sed 's/:.*//g'`
IN=($INIT)
EN=($END)
for P in 0 1 2
	do
		I=${IN[$P]} 
		E=${EN[$P]} 
		cat $FILE | sed -n "${I},${E}p" | grep -v "Mutations" | grep -v "Genomes" > ${FILE}_pop${P}
		head ${FILE}_pop${P}
	done

# P 0 is p1 outgroup
# P 1 is p2 wild
# P 2 is p3 domesticated