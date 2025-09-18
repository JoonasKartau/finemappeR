README
================

## Installation

The R package can be installed using the following command in R.

``` r
# install.packages("remotes")  # Uncomment if remotes is not installed
remotes::install_github("JoonasKartau/finemapmiss", auth_token = NULL)
```

You may still need a GitHub PAT token. (Steps copied from https://github.com/orgs/community/discussions/140956)

## Setting up a GitHub Personal Access Token (PAT) in R

1. Go to `Github/settings/tokens`.  
2. Click **"Generate new token." (classic)**.  
3. Select the appropriate scopes for your use case. For installing packages, the following scopes are usually enough:  
   - `repo`  
   - `read:packages`  
4. Click **"Generate token."**  
5. Copy the token, as you'll need it soon.

### Set up the PAT in R

1. Open R and run the following command to set the new token:

``` r
gitcreds::gitcreds_set()
```

## Example run

Once installed, an example of meta-analysis fine-mapping can be accessed as a vignette:

```r
vignette("FINEMAP-miss", package = "finemapmiss")
```

