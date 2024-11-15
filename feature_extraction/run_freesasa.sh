pdb_sort $1 | pdb_tidy | pdb_reatom | pdb_delchain -C | pdb_sort | pdb_tidy | pdb_reatom | grep -E 'ATOM  |TER  |END  ' > "${2}${3}_receptor.pdb"
touch "${2}${3}_rsa_sasa.csv"
peptide_length=$(echo -n "$4" | wc -c)
for ((i = 1; i <= ${peptide_length}; ++i)); do
	pdb_sort $1 | pdb_tidy | pdb_reatom | pdb_delchain -A | pdb_sort | pdb_tidy | pdb_reatom | grep -E 'ATOM  |TER  |END  ' | pdb_selres -${i} > "${2}${3}_residue_${i}.pdb"
	pdb_merge "${2}${3}_receptor.pdb" "${2}${3}_residue_${i}.pdb" | pdb_sort | pdb_tidy | pdb_reatom > "${2}${3}_complex.pdb"
	freesasa "${2}${3}_complex.pdb" -n 20 --format=rsa --radii=naccess >> "${2}${3}_freesasa.log"
done
cat "${2}${3}_freesasa.log" | grep "C  " | grep "RES" | awk '{print $2, $4, $5, $6, $7, $8}' | tr ' ' ',' | awk -F, -v pep="$index", -v OFS="," '{$7=pep}1' | awk -F',' '{print $7, $1, $2, $3, $4, $5, $6}' | tr ' ' ',' >> "${2}${3}_rsa_sasa.csv"