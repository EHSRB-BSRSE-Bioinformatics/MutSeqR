# # API calls
# 
# API method: /file-xxxx/download

# Specification
# Generates a "download URL" for downloading the contents of this file object. The download URL may refer to a different endpoint than the DNAnexus API server, and accepts HTTP GET requests.
# Requests to the download URL must be initiated within the number of seconds specified in the "duration" input parameter (starting from the time this call is made, according to the server), after which the URL will expire. GET requests MUST include any headers specified in the API server's response to /file-xxxx/download (see below). The download URL also honors the "Range" HTTP request headers, allowing clients to download only a particular byte range of the file.
# The download URL has the following support for CORS:
# If a GET request to download URL includes the "Origin" header, its contents will be propagated into the "Access-Control-Allow-Origin" header of the response.
# Preflight requests (OPTIONS requests to a part upload URL, with appropriate extra headers as defined in the CORS draft) will be accepted if the value of the "Access-Control-Request-Method" header is "GET". The values of "Origin" and "Access-Control-Request-Headers" (if any) of the request, will be propagated to "Access-Control-Allow-Origin" and "Access-Control-Allow-Headers" respectively in the response. The "Access-Control-Max-Age" of the response is set to 1 hour.
# Successful calls to the download URL will return the HTTP response code 200, and will include a "Content-Type" header, set to whatever Internet Media Type was specified when the file object was created, and a "Content-Disposition: attachment" header that may also include a filename, if requested (see below). The request may include the query string "?inline" to override the Content-Disposition header. Unsuccessful requests will return an HTTP error response code (and in that case there are no guarantees about the response body, as the download URL does not necessarily conform to the general API rules regarding error messages).

# Inputs
# duration int (optional, default is equivalent to 24 hours) Number of seconds (starting from the time this call is made, according to the server) during which the generated URL will be valid
# filename string (optional) The desired filename of the downloaded file, to be affixed to the returned URL.  If provided, this filename will be encoded as a URI component and affixed to the download URL, whose resource portion will end in e.g. '/filename', to ease downloads through web browsers and utilities such as wget.
# project string (optional) ID of a project containing the file, with which the download URL will be associated. Requests to the download URL will succeed only so long as 1. the file still resides in this project and 2. the user who generated the URL still has at least VIEW permission to this project. If this value is not provided, the URL will work so long as the file resides in any project to which the user who generated the URL has at least VIEW permission. This field must be provided to get the download URL for a watermarked file when invoked outside the context of a DNAnexus job.
# preauthenticated boolean (optional, default false) Whether to generate a "preauthenticated" download URL, which embeds any necessary authentication information in the URL itself, rather than requiring separate request headers
# Security note: URLs generated in this way intrinsically provide access to the file data to anyone in possession. Therefore, they should not be unnecessarily stored, logged, printed to console, etc. in production applications.
# For security reasons, preauthenticated URLs should be project specific.
# stickyIP boolean (optional if preauthenticated is true; required to be false otherwise; default false) Whether HTTP GET requests to the preauthenticated download URL should be restricted to a single origin IP address. If stickyIP and preauthenticated are true, then the first HTTP GET request to the preauthenticated download URL will dictate the IP address from which all subsequent requests must originate.

# Outputs
# url string An absolute URL to which HTTP GET requests can be made to download the file
# headers mapping HTTP headers which MUST be supplied with any GET request to url
# key Header field name
# value string Header value
# Security note: the headers may include authentication tokens, and therefore should not be stored, logged, printed to console, etc. in production applications.
# Note: if a preauthenticated URL was requested, then no keys will be present.

# Errors
# ResourceNotFound
# project is specified but the file object is not in the specified project
# PermissionDenied
# VIEW access required to some project that contains the file object
# If project is specified, VIEW access is required to that project
# InvalidInput
# duration (if provided) is not a positive integer
# InvalidState
# The file object is not in the "closed" state

# File downloads in web applications
# To generate non-preauthenticated file download URLs, web applications (running inside web browsers) should make /file-xxxx/download requests to the separate endpoint https://dl.dnanex.us instead of https://api.dnanexus.com. Browser requests to non-preauthenticated file download URLs are authenticated by means of a URL-specific cookie, set by the API server's response to the /file-xxxx/download route on this separate endpoint.
# Non-browser-based applications implementing the above specification, or web applications only needing preauthenticated download URLs, may call /file-xxxx/download on https://api.dnanexus.com as usual.

library("httr")
token <- Sys.getenv("DNANexus_API_key")
GET("https://api.dnanexus.com/system/whoami", query = list("Bearer" = token))



# https://auth.dnanexus.com/oauth2/authorize

# https://dnanexus.com/system/whoami

# Example: Authorization: Bearer 7Fjfp0ZBr1KtDRbnfVdmIw

