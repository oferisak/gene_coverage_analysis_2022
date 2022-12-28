library(usethis)

# Step 1. setup git in project
setup_git<-function(){
  ## now, working inside "the package", initialize git repository
  use_git()
}

# Step 2. get git token
create_git_token<-function(){
  ## create github token
  create_github_token()
  # set it up
  gitcreds::gitcreds_set()
}

# Step 3. create git repository
create_git_repo<-function(){
  ## create github repository and configure as git remote
  use_github()
}

# After these steps are taken, you can push/pull from github