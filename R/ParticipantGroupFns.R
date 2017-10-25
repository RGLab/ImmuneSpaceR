###############################################################################
###              PARTICIPANT FILTERING METHODS                              ###
###############################################################################

.ISCon$methods(
  listParticipantGroups = function(){
    "returns a dataframe with all saved Participant Groups on ImmuneSpace.\n"
    
    if(config$labkey.url.path != "/Studies/"){
      stop("labkey.url.path must be /Studies/. Use CreateConnection with all studies.")
    }
    
    user <- .validateUser(con = .self)
    
    # Get summary data
    pgrpSql <- paste0("SELECT ",
                      "ParticipantGroup.RowId AS GroupId, ",
                      "ParticipantGroup.Label AS GroupName, ",
                      "COUNT(*) AS Subjects ",
                      "FROM ",
                      "ParticipantGroupMap ",
                      "FULL OUTER JOIN ParticipantGroup ",
                      "ON ParticipantGroupMap.GroupId=ParticipantGroup.RowId ",
                      "WHERE ",
                      "ParticipantGroup.CategoryId IN (",
                      "SELECT rowId ",
                      "FROM ParticipantCategory ",
                      "WHERE OwnerId = (",
                      "SELECT Users.UserId ",
                      "FROM core.Users ",
                      "WHERE Users.Email = ",
                      sprintf("'%s'", user), " )", # need single quotes for string
                      ") ",
                      "GROUP BY ",
                      "ParticipantGroup.RowId, ",
                      "ParticipantGroup.Label ")
    
    result <- labkey.executeSql(config$labkey.url.base,
                                config$labkey.url.path,
                                schemaName = "study",
                                sql = pgrpSql,
                                colNameOpt = "fieldname",
                                showHidden = TRUE)
    
    if( nrow(result) == 0 ){
      warning(paste0("No participant groups found for user email: ", user))
    }
    
    return(result)
  }
)

.ISCon$methods(
  getParticipantData = function(group, dataType, original_view = FALSE, ...){
    "returns a dataframe with ImmuneSpace data subset by groupId.\n
    group: Use con$listParticipantGroups() to find Participant groupId or groupName.\n
    dataType: Use con$listDatasets('datasets') to see possible dataType inputs.\n"
    
    if(config$labkey.url.path != "/Studies/"){
      stop("labkey.url.path must be /Studies/. Use CreateConnection with all studies.")
    }
    
    if(!(dataType %in% .self$available_datasets$Name)){
      stop("DataType must be in ", paste(.self$available_datasets$Name, collapse = ", "))
    }
    
    # Checking to ensure user is owner of group on correct machine
    # Note that groupId is used for actually sql statement later regardless of user input
    user <- .validateUser(con = .self)
    
    # Must rerun this to ensure valid groups are only for that user and are not changed
    # within cache.
    validGrps <- .self$listParticipantGroups()
    
    if(typeof(group) == "double"){
      col <- "GroupId"
      groupId <- group
    }else if( typeof(group) == "character"){
      col <- "GroupName"
      groupId <- validGrps$GroupId[ validGrps$GroupName == group ]
    }else{
      stop("Group ID or Name not interpretable as double or character. Please reset")
    }
    
    if( !( group %in% validGrps[[col]] ) ){
      stop(paste0("group ", group,
                  " is not in the set of ", col,
                  " created by ", user,
                  " on ", .self$config$labkey.url.base))
    }
    
    dt <- dataType # for brevity
    
    # Get assay data + participantGroup + demographic data
    # Handle special cases
    if( dt == "demographics"){
      varSelect <- "cohort_membership.name AS Cohort "
      varJoin <- paste0(" JOIN cohort_membership ",
                        "ON ",
                        dt, ".ParticipantId = cohort_membership.ParticipantId ")
    }else{
      varSelect <- paste0("demo.age_reported, ",
                          "demo.race, ",
                          "demo.gender, ")
      varJoin <- paste0(" JOIN immport.arm_or_cohort ",
                        "ON ", 
                        dt, ".arm_accession = arm_or_cohort.arm_accession ")
      if( dt %in% c("gene_expression_files", "cohort_membership", "fcs_sample_files" ) ){
        varSelect <- paste0( varSelect,
                             "demo.age_unit, ",
                             "demo.age_event, ",
                             "demo.ethnicity, ",
                             "demo.species, ",
                             "demo.phenotype ")
      }else{
        varSelect <- paste0( varSelect, "arm_or_cohort.name AS arm_name ")
      }
    } 
    
    sqlAssay <- paste0("SELECT ",
                       dt, ".*, ",
                       "pmap.GroupId AS groupId, ",
                       varSelect,
                       "FROM ",
                       dt,
                       " JOIN ParticipantGroupMap AS pmap ",
                       "ON ", dt, ".ParticipantId = pmap.ParticipantId",
                       " JOIN Demographics AS demo ",
                       "ON ", dt, ".ParticipantId = demo.ParticipantId",
                       varJoin,
                       " WHERE pmap.GroupId = ", as.character(groupId))
    
    assayData <- labkey.executeSql(baseUrl = .self$config$labkey.url.base,
                                   folderPath = .self$config$labkey.url.path,
                                   schemaName = "study",
                                   sql = sqlAssay,
                                   colNameOpt = "fieldname",
                                   ...) # allow for params to be passed, e.g. maxRows
    
    # Want to match getDataset() results in terms of colnames / order
    defaultCols <- colnames(.self$getDataset(x = dt,
                                             original_view = original_view,
                                             maxRows = 1)) 
    
    # Some names from assayData do not match default cols and need to changed manually.
    colnames(assayData)[ grep("ParticipantId", colnames(assayData)) ] <- "participant_id"
    
    if( original_view == FALSE){
      changeCol <- if(dt == "demographics"){ 
        "Cohort"
      }else if(dt == "cohort_membership"){ 
        "name"
      }else{ 
        "arm_name"
      }
      
      colnames(assayData)[ grep(changeCol, colnames(assayData)) ] <- "cohort"
    }else{
      if( dt %in% c("gene_expression_files", "fcs_control_files") ){
        colnames(assayData)[ grep("arm_name", colnames(assayData)) ] <- "cohort"
      }
    }
    
    filtData <- assayData[ , colnames(assayData) %in% defaultCols ]
    filtData <- filtData[ , order( match(colnames(filtData), defaultCols)) ]
    
    return(filtData)
  }
)