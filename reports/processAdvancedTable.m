function processAdvancedTable(doc, reportData, holeId)
   %This function adds a table to our report, which is solely constructed
   %with DOM API commands
   import mlreportgen.dom.*
   
   %We add a new Level2 Heading before this table
   p = Paragraph('Complex Table Example', 'AR_Heading2');
   append(doc, p);
   %And some description
   p = Paragraph(['This is an example of an advanced table, which is ' ...
      'completely created with the DOM API.'], 'AR_Normal');
   append(doc, p);
   
   tableData = reportData.(holeId);
   numCol = size(tableData,2);
   numRow = size(tableData,1);
   %The FormalTable class allows us to build a complete table from scratch.
   %We need to provide the number of columns in the constructor
   table = FormalTable( numCol );
   %The table shall span the page width, so we set the attribute to 100%
   table.Width = '100%';

   %Now we need to create a TableRow object for each row we want to add.
   %This loop is for adding a table header which shall be displayed in bold
   row = TableRow();
   for nc=1:numCol
      %We create a Text for each colum in the header and make it bold
      textObj = Text(sprintf('Row %d', nc));
      textObj.Bold = true;      
      te = TableEntry( textObj );
      %Then a TableEntry is added to the TableRow for each column
      append(row, te );
   end
   %This row is appended to the table-header
   append(table.Header, row);
   
   %This loop fills the table body
   for nr=1:numRow
      row = TableRow();
      for nc=1:numCol
         te = TableEntry( tableData{nr,nc} );
         %The BackgroundColor of each TableEntry shall be filled with a
         %random color.
         bgColor = sprintf('#%x', randi(16777215));
         te.Style  = { BackgroundColor(bgColor) };
         te.VAlign = 'middle';
         append(row, te );
      end
      append(table.Body, row);
   end   
   append(doc, table);   
end
