from __future__ import annotations

import pandas as pd
from sqlalchemy.engine import Engine

from ..sql.models import Base
from ..logger import log


def select_from_table(
    engine: Engine, tablename: str, column: str, value: str
) -> pd.DataFrame | None:
    try:
        df = pd.read_sql_table(tablename, engine)
        result = df[df[column] == value]
    except (LookupError, ArithmeticError, TypeError) as e:
        log.debug("unable to fetch from database table", exc_info=e)
        result = None

    return result


def create_alert_tables(engine: Engine):
    Base.metadata.create_all(engine)
