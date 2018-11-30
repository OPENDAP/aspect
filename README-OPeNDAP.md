URL File Reader extension
====================

About
-----

The utilities.cc file has been modified to read data from the OPeNDAP servers. 
If a URL is passed in as the filename, the read_and_distribute_file_content() will read the data 
from the provided URL and treat it as it would a local file.

URL File Format
---------------

The url reader requires data to be read in a certain way, but there are a few different ways of
storing that data for Aspect to get. We used CSV files to store and read our data, but other
methods, such as NetCDF, can also be used. The data must be formated properly so that our
url reader can obtain the necesary information (column header, POINTS, data type, and data).

The function can read any number of columns of data as well as any number of rows. The "# POINTS:" 
values are represented as attributes in the CSV file. The attributes can be set using a .das file of 
the same name as the .csv file. Both the .csv and the .csv.das file must be included 
(i.e. **urlFile.csv** and **urlFile.csv.das**).

CSV files must include the column title followed by <String>, this is because Aspect takes in 
values from local files as strings (and url files are designed to be read in to look like local 
files). The CSV must then be organized so that there are the same number of values per row as 
there is number of columns. So 4 columns would mean that there has to be 4 values and then a 
new line.
	

Example CSV
-----------

	"latitude<String>", "longitude<String>", "depth<String>"
	0,1,2
	5,4,3

Example DAS
-----------

	Attributes {
	    depth {
	        int32 points 2;
	    }
	}
