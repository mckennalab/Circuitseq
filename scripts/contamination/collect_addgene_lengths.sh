date=$(date '+%Y-%m-%d')

#find addgene_simulated_data/ -type f -name "*.fa" -print0 | tail -n +2 | wc --files0-from=-| sed 's/addgene_simulated_data\/[[:digit:]]\+_test\///g' | sed 's/\.fa//g' | grep -v otal > $(date '+%Y-%m-%d')_collected_addgene_lengths.txt

find addgene_simulated_data/ -type f -name "*.fa" | xargs -I {} sh -c 'tail -n +2 {} | wc -cl' > $(date '+%Y-%m-%d')_collected_addgene_lengths.txt
find addgene_simulated_data/ -type f -name "*.fa" | sed 's/addgene_simulated_data\/[[:digit:]]\+_test\///g' | sed 's/\.fa//g' > $(date '+%Y-%m-%d')_collected_addgene_names.txt
