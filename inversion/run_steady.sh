mkdir $1/output/
mkdir $1/steady/
seq 100 | parallel -j44 "python inversion/kf/run_steady.py $1 {}  >> $1output/steady_{%}.txt"
