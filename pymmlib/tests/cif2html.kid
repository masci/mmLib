<?xml version='1.0' encoding='utf-8'?>

<html xmlns="http://www.w3.org/1999/xhtml" xmlns:py="http://purl.org/kid/ns#">
  <head>
    <title>${cif.struct.entry_id}</title>
  </head>
  <body>
    <center><h1>
      <a href="http://www.rcsb.org/pdb/explore.do?structureId=${cif.struct.entry_id}">
      ${cif.struct.entry_id}</a>
    </h1></center>
    <p style="font-size:large">${cif.struct.title}</p>
    <p style="font-size:large">${cif.struct.pdbx_descriptor}</p>

    <table class="summary">
      <tr><td>${cif.exptl.method}</td></tr>
      <tr><td>${cif.refine.ls_d_res_high}</td></tr>
      <tr><td>${cif.refine.pdbx_isotropic_thermal_model}</td></tr>
      <tr><td>"${cif.symmetry["space_group_name_H-M"].replace(" ","")}"</td></tr>
      <tr><td>${cif.citation.title}<br/>${cif.citation.year}</td></tr>
      <tr><td>${cif.struct_keywords.pdbx_keywords}</td></tr>
      <tr><td>${cif.struct_keywords.text}</td></tr>
    </table>

    <table class="authors">
      <tr py:for="ca in cif.citation_author">
        <td>${ca.name}</td>
      </tr>
    </table>

    <table class="dates">
      <tr py:for="rev in cif.database_PDB_rev">
        <td>${rev.date}</td>
      </tr>
    </table>

    <table class="software">
      <tr py:for="software in cif.software">
        <td>${software.name}</td>
      </tr>
    </table>

  </body>
</html>
