<?xml version='1.0' encoding='utf-8'?>

<html xmlns="http://www.w3.org/1999/xhtml" xmlns:py="http://purl.org/kid/ns#">

  <head>
    <title>TLSMD Hinge Predictions for Structure ID ${tlsmd.struct_id}</title>
    <link rel="stylesheet" href="http://www.drizzle.com/~jpaint/molmovdb.css" type="text/css" media="screen"/>
  </head>

  <body>
    <div id="page">
      <center><h1>TLSMD Hinge Predictions for Structure ID ${tlsmd.struct_id}</h1></center>

      <center>
      <div py:for="chain in tlsmd.iter_chains()">
	<table id="partitiontable" width="75%">
	  <div py:for="title in chain.tbl.iter_column_titles()">
            <th>${title}</th>
          </div>

          <tr py:for="i, row in enumerate(chain.tbl.iter_rows())" class="${not int(row[2]) and 'selectedrow' or  i%2 and 'oddrow' or 'evenrow'}">
            <div py:for="item in row">
              <td>${item}</td>
            </div>
          </tr>
	</table>
      </div>
      </center>

    </div>
  </body>
</html>
