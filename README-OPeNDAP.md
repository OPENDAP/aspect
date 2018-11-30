URL File Reader extension
====================

About
-----

Aspect has been modified so that it can read data from OPeNDAP data
servers in addition to local files. OPeNDAP is a RESTful web protocol
used to access remote data. Several groups provide data servers that
support the protocol. For more information, see www.opendap.org. This
work was done as part of the BALTO NSF EarthCube grant (# ...).

For Aspect, if a URL is used in place of a file name in a prm file, the function
_read\_and\_distribute\_file\_content()_ will read the data from the
provided URL and treat it as it would a local file. The changes to the
Aspect software are limited to the utilities.cc source file and to the
cmake scripts. See INSTALL-OPeNDAP for more information about
building Aspect so that it includes support for OPeNDAP.

Serving Remote Data
---------------

The OPeNDAP URL reader requires data to be presented to Aspect a
certain way, although there are a number of different ways of storing
those data at the server. In essence, the server provides Aspect with
a uniform interface to use for reading data. How the data are
organized is primarily an issue for the data provider.

To test this new functionality, we used Lithospheric Thickness data
stored in CSV files. The Hyrax OPeNDAP server can easily serve this
kind of data and the sample data were already stored in those kinds of
files. We added some additional information to provide column headers
and the number of data points. The "# POINTS:" values are represented
as dataset attributes stored in an ancillary file. The attributes can
be set using a file the extension _.das_ and the same base name as the
_.csv_ data file file. Both the .csv and the .csv.das file must be in
the same directory for this to work.(i.e. **urlFile.csv** and
**urlFile.csv.das**). You can see these data at
http://test.opendap.org/opendap/BALTO/lithospheric_thickness.csv.html

While the sample data used has four columns and about 240 row, the
_read\_and\_distribute\_file\_content()_ function can read any number
of columns of data as well as any number of rows.

Note that these data could have beeen stored in an SQL database, a
netCDF file or other data storage technology so long as an OPeNDAP
server could read it and provide the information in the way the Aspect
software expects.

Example CSV File
-----------

	"longitude<String>", "colatitude<String>", "litho<String>", "crust<String>"

    "0.261799", "1.22173", "61630", "43038"
    "0.266163", "1.22173", "61490", "41447"
	...

Example DAS File
-----------

	Attributes {
	    depth {
	        int32 points 2;
	    }
	}
