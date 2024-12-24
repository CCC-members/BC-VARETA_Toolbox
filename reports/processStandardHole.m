function processStandardHole(doc, reportData, holeId)
   import mlreportgen.dom.*
   
   %This function fills all holes that do not require a special handling.
   %The try/catch block performs an error handling in the case that there
   %is a mismatch between the hole-names in the template and the field-
   %names of the reportData-structure
   try
      data = reportData.(holeId);
   catch
      warning('Undefined Hole-Id: %s', holeId);
      data = 'undefined';
   end

   %Create a new text object and add it to the report
   t = Text(data);   
   append(doc, t);
end

