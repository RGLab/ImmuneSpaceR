# Set curlOptions for download ---------------------------------
if (!any(file.exists("~/.netrc", "~/_netrc"))) {
  ISR_login <- Sys.getenv("ISR_login")
  ISR_pwd <- Sys.getenv("ISR_pwd")
  ISR_machine <- ifelse(Sys.getenv("ISR_machine") == "",
                        "www.immunespace.org",
                        Sys.getenv("ISR_machine"))
  if (ISR_login != ""  &  ISR_pwd != "") {
    netrc_file <- tempfile("ImmuneSpaceR_tmp_netrc")
    netrc_string <- paste("machine", ISR_machine, 
                          "login", ISR_login, 
                          "password", ISR_pwd)
    write(x = netrc_string, file = netrc_file)
    labkey.netrc.file <- netrc_file
  }
}