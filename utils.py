from sqlalchemy import create_engine

def connectDB():
    db_connection_url = "postgresql://postgres:wave*worm@88.99.137.51:5432/maceio_tests"
    engine = create_engine(db_connection_url)
    return engine