# Set curlOptions for download ---------------------------------
ISR_login <- Sys.getenv("ISR_login")
ISR_pwd <- Sys.getenv("ISR_pwd")
if(ISR_login != ""  &  ISR_pwd != ""){
  netrc_file <- tempfile("ImmuneSpaceR_tmp_netrc")
  netrc_string <- paste("machine www.immunespace.org login", ISR_login, "password", ISR_pwd)
  write(x = netrc_string, file = netrc_file)
  labkey.netrc.file <- netrc_file
}