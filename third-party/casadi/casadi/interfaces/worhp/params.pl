while (<>) {
  $description = $1 if /\* (.*) \*/;
  $type = "OT_BOOL" if /^\s+bool/;
  $type = "OT_DOUBLE" if /^\s+double/;
  $type = "OT_INT" if /^\s+int/;

  if (/^\s+(\w+)\s+(\w+)/ and $2 ne"Ares") {
    print qq|addOption("$2",$type,worhp_p.$2,"$description");\n|;
    #print qq|if (hasSetOption("$2")) worhp_p.$2 = option("$2");\n|;
  }

}
