for i in *.dna ; do filename=${i%.dna}
		    python read_dna.py $i $filename'_ref.fasta'
		    python clean_ref.py $filename
		    rm $filename'_ref.fasta'
		    done
		    
