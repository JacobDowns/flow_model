mkdir $1/output/
mkdir $1/sigmas/
seq 2000 | parallel -j44 "python inversion/ut/run_sigma.py $1 {} >> $1output/{%}.txt"
python3 inversion/ut/write_sigmas.py $1/
