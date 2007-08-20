
function ToggleDivVisibility(control_id, target_id, show_val, hide_val) {  
  
  var ctrl_element = document.getElementById(control_id); 
  var target_element = document.getElementById(target_id);

  if (target_element.style.display != "none") {
    target_element.style.display = "none";
    ctrl_element.firstChild.nodeValue = show_val;
  } else {
    target_element.style.display = "inline";
    ctrl_element.firstChild.nodeValue = hide_val;
  }
}
