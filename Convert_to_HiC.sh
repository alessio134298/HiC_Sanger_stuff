# Takes in input two arguments: $1 that is the matrix, $2 s the chromssizes file and $3 that can be "h5" or "mcool"

# bash Convert_to_hic.sh matrix.h5 chromsizes.txt h5

# Note: in juicer-pre command (the last one that performs the final conversion to .hic) leave the max amount of memory usage to 128 GB,
# with less than this it crashes due to insufficient memory

module load HGI/softpack/groups/escramble/eSCRAMBLE/8

echo "Input is: $1"

FOLDER=$(dirname $1)

if [ "$3" == "h5" ]
then
	MATRIX=$(/usr/bin/basename "$1" | /usr/bin/awk -F ".h5" '{print $1}')

	echo "Output directory will be: ${FOLDER}"
	echo "Name of the file is: ${MATRIX}"

	# .h5 to mcool (useful to have the mcool format)

	echo "Conversion to mcool"

	hicConvertFormat \
		-m  ${FOLDER}/${MATRIX}.h5 \
		--inputFormat h5 \
		--outputFormat mcool \
		-r 1000 2500 5000 10000 20000 50000 100000 250000 500000 1000000 \
		-o ${FOLDER}/${MATRIX}.mcool

	# .h5 to .cool

	echo "Conversion to cool"

	hicConvertFormat \
		-m ${FOLDER}/${MATRIX}.h5 \
		--inputFormat h5 \
		--outputFormat cool \
		-o ${FOLDER}/${MATRIX}.cool \
		--resolutions 1000 &&

	# convert to ginteractions

	echo "Conversion to ginteractions"

	hicConvertFormat \
		-m ${FOLDER}/${MATRIX}.cool \
		-o ${FOLDER}/${MATRIX}.ginteractions \
		--inputFormat cool --outputFormat ginteractions &&

	# format resulting file
	bash /lustre/scratch126/gengen/teams/parts/ab77/scripts/awk.sh \
	${FOLDER}/${MATRIX}.ginteractions.tsv \
	${FOLDER}/${MATRIX}.ginteractions.short.tsv &&

	# sort the file
	sort -k2,2d -k6,6d \
	${FOLDER}/${MATRIX}.ginteractions.short.tsv \
	> ${FOLDER}/${MATRIX}.ginteractions.sorted.tsv &&

	mkdir -p ${FOLDER}/tmp &&

	# convert
	
	echo "Conversion to hic"

	java -Xms4g -Xmx128g -jar /lustre/scratch126/gengen/teams/parts/ab77/scripts/packages/juicer_tools_1.22.01.jar pre -j 12 \
		-t ${FOLDER}/tmp \
		-r 1000,2500,5000,10000,20000,50000,100000,250000,500000,1000000 \
		${FOLDER}/${MATRIX}.ginteractions.sorted.tsv \
		${FOLDER}/${MATRIX}.hic \
		$2

	rm ${FOLDER}/${MATRIX}.ginteractions.tsv ${FOLDER}/${MATRIX}.ginteractions.short.tsv ${FOLDER}/${MATRIX}.ginteractions.sorted.tsv

elif [ "$3" == "mcool" ]
then
	MATRIX=$(/usr/bin/basename "$1" | /usr/bin/awk -F ".mcool" '{print $1}')

	echo "Output directory will be: ${FOLDER}"
	echo "Name of the file is: ${MATRIX}"

	# .mcool to .cool
	cooler cp ${FOLDER}/${MATRIX}.mcool::resolutions/1000 ${FOLDER}/${MATRIX}.cool 

	# convert to ginteractions
	hicConvertFormat \
		-m ${FOLDER}/${MATRIX}.cool \
		-o ${FOLDER}/${MATRIX}.ginteractions \
		--inputFormat cool --outputFormat ginteractions &&

	# format resulting file
	bash /lustre/scratch126/gengen/teams/parts/ab77/scripts/awk.sh ${FOLDER}/${MATRIX}.ginteractions.tsv ${FOLDER}/${MATRIX}.ginteractions.short.tsv &&

	# sort the file
	sort -k2,2d -k6,6d \
	${FOLDER}/${MATRIX}.ginteractions.short.tsv \
	> ${FOLDER}/${MATRIX}.ginteractions.sorted.tsv &&

	mkdir -p ${FOLDER}/tmp &&

	# convert
	java -Xms4g -Xmx128g -jar /lustre/scratch126/gengen/teams/parts/ab77/scripts/packages/juicer_tools_1.22.01.jar pre -j 12 \
		-t ${FOLDER}/tmp \
		-r 1000,2500,5000,10000,20000,50000,100000,250000,500000,1000000 \
		${FOLDER}/${MATRIX}.ginteractions.sorted.tsv \
		${FOLDER}/${MATRIX}.hic \
		$2

	rm ${FOLDER}/${MATRIX}.ginteractions.tsv ${FOLDER}/${MATRIX}.ginteractions.short.tsv ${FOLDER}/${MATRIX}.ginteractions.sorted.tsv	
fi