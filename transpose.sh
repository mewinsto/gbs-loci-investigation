##Transpose 
#lines = wc outfile | cut -f 2 -d " "
#bytes = wc outfile | cut -f 4 -d " "
#width = lines/bytes
#in this case width = 48
COUNTER=1
while [   $COUNTER -lt 48 ]; do
    cut -b $COUNTER outfile >> hohoho.txt
    let COUNTER=COUNTER+1
done
