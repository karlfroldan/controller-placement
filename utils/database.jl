module Utils_Database

# using SQLite

import SQLite

export Database

const DEFAULT_DB_STMTS :: String =
    """
    CREATE TABLE strategy (
        id INTEGER PRIMARY KEY,
        
    );
    # CREATE TABLE iterations    
    """

mutable struct Database
    db :: SQLite.DB
end

function load_database(dbname :: String) Database
    # db :: SQLite.DB = undef
    # SQLite.load!(
    db = SQLite.DB(dbname)
    Database(db)
end

end
