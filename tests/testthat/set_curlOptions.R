# Set curlOptions for download ---------------------------------
if (!any(file.exists("~/.netrc", "~/_netrc"))) {
  labkey.netrc.file <- ImmuneSpaceR:::get_env_netrc()
  labkey.url.base <- ImmuneSpaceR:::get_env_url()
}
