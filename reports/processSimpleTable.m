function processSimpleTable(doc, reportData, holeId)
   %This function adds a simple table, which is based on the Word style
   %"_AR_Table" to our report
   import mlreportgen.dom.*
   
   %We add a new Level2 Heading before this table
   p = Paragraph('Simple Table Example', 'AR_Heading2');
   append(doc, p);
   %And some description
   p = Paragraph(['This is an example of a simple table. The table ' ...
      'is based upon the Word table style "AR_Table".'], 'AR_Normal');
   append(doc, p);
   
   tableData = reportData.(holeId);
   %The Table class constructor accepts a cell array which contains the
   %data we want to display. The complete table is automatically build for
   %us, based on the input data and the table style AR_Table, which is
   %defined in the Word template
   table = Table( tableData, 'AR_Table');
   append(doc, table);   
end
