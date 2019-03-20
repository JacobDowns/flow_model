mkdir $1/output/
seq 700 | parallel -j46 "python inversion/ut/run_sigma.py $1 {} >> $1output/{%}.txt"
python3 inversion/ut/write_sigmas.py $1/
