
# Extract substrings to check for primer matching

mismatch=1

while getopts i:p:n:s: option
do 
    case "${option}"
        in
        i)input_file=${OPTARG};;
        s)strand=${OPTARG};;
    esac
done

if [ $strand == 'F' ]; then
    zcat $input_file | sed -n '2~4p' | awk '{print substr($0,1,50)}' > ./substrings.txt
elif [ $strand == 'R' ]; then
    zcat $input_file | sed -n '2~4p' | awk '{print substr($0,length($0)-49,50)}' > ./substrings.txt
fi