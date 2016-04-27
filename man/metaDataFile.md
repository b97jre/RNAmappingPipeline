#Meta data file 

## Initial setup of metadata file

|SampleName     |Name  |Prep_date|Seq_Date|Seq_type|Sep1|Sep2|Sep3|forward|reverse|
|--|--|--|--|--|--|--|--|--|--|--|--|
|SampleName1     |SN1 |1|1|Ilumina|1|1|1|/path/to/forwardFile1 |/path/to/reverseFile1|
|SampleName3     |SN2 |1|2|Solid|2|2|1|/path/to/forwardFile2 |/path/to/reverseFile2|
|SampleName3     |SN2 |2|3|PacBio|3|1|2|/path/to/forwardFile3 |/path/to/reverseFile3|

## Examples of filled in meta data file



### Example with three samples from two different tissues from three persons with the same mother but different fathers

|SampleName|Name|Prep_date|Seq_Date|Seq_type|Gender|Individual|Mother|Father|Operation_date|Tissue|weight|forward|reverse|	
|----------|----|---------|--------|--------|------|----------|------|------|--------------|------|------|------|------|	
|SampleName1 |SN1 |1       |1      |Ilumina |Male |1| 1|1|1|Tissue1 |86.7 |/path/to/forwardFile1 |/path/to/reverseFile1|
|SampleName2 |SN2 |1       |1      |Ilumina |Male |2| 1|2|2|Tissue1 |88.7 |/path/to/forwardFile2 |/path/to/reverseFile2|
|SampleName3 |SN3 |1       |1      |Ilumina |Female |3| 1|1|3|Tissue1 |74.7 |/path/to/forwardFile3 |/path/to/reverseFile3|
|SampleName4 |SN1 |1       |1      |Ilumina |Male |1| 1|1|1|Tissue2 |86.7 |/path/to/forwardFile4 |/path/to/reverseFile4|
|SampleName5 |SN2 |1       |1      |Ilumina |Male |2| 1|2|2|Tissue2 |88.7 |/path/to/forwardFile5 |/path/to/reverseFile5|
|SampleName6 |SN3 |1       |1      |Ilumina |Female |3| 1|1|3|Tissue2 |74.7 |/path/to/forwardFile6 |/path/to/reverseFile6|




