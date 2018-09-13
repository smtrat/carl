if [ `uname` = 'Linux' ]; then
	find . -type f | xargs sed -i "s/auto_ptr/unique_ptr/g"
elif [ `uname` = 'Darwin' ]; then
	export LC_ALL=C
	find . -type f | xargs sed -i '' "s/auto_ptr/unique_ptr/g"
else
	echo "Could not patch CoCoA, unknown system type `uname`"
fi