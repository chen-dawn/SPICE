# How to prepare the data file

first download the things that are in this [Dropbox folder](https://www.dropbox.com/scl/fo/qqazqiqemam34uqt1abxy/h?rlkey=48v58t32kaxjgkyewhth1ph7m&e=1&dl=0).

Then make the dataset and data file by:

```
python create_datafile.py train all
python create_datafile.py test 0    

python create_dataset.py train all  
python create_dataset.py test 0
```

