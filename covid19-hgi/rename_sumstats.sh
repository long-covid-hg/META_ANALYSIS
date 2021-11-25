for i in ls Analysis_*.gz
do

ana="${i:9:2}"

Nca=$(zcat ${i} | head -n 2 | awk 'FNR==2 {print $21}')
Nco=$(zcat ${i} | head -n 2 | awk 'FNR==2 {print $22}')

if [[ ${i} == *"_M.gz"* ]]
then
    sex=M
    age=ALL
elif [[ ${i} == *"_F.gz"* ]]
then
    sex=F
    age=ALL
elif [[ ${i} == *"_GT_60.gz"* ]]
then
    sex=ALL
    age=GT_60
elif [[ ${i} == *"_LE_60.gz"* ]]
then
    sex=ALL
    age=LE_60
else
    sex=ALL
    age=ALL
fi

cp ${i} GENCOVID.Pigazzini.ANA_${ana}_V2.2.${age}.${sex}.EUR.${Nca}.${Nco}.SAIGE.20200927.txt.gz
done