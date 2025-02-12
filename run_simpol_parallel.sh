output_file="simpol_commands.txt"

> "$output_file"

for k in 20 30 40 50 60 80 100; do
    for kSd in 0 5 10; do
        for a in 1 1.5 2 3.5 5; do
            for b in 0.5 1 1.5 2; do
                for z in 1500 1750 2000 2250 2500; do
                    for zetaSd in 500 1000 1500; do
                        for t in 0 15 30 45 60 75 90 105 120; do
                            output_dir="${k}k${kSd}kSd${a}a${b}b${z}z${zetaSd}zetaSd${t}t_results"
                            echo "./Release/simPol -k $k --kSd $kSd --kMin 15 --kMax 200 --geneLen 3000 -a $a -b $b -z $z --zetaSd $zetaSd --zetaMax 2500 --zetaMin 1500 -n 100 -s 33 --addSpace 17 -t $t -d $output_dir" >> "$output_file"
                        done
                    done
                done
            done
        done
    done
done

./ParaFly -c simpol_commands.txt -CPU 10