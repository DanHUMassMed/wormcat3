from enum import Enum


class ErrorCode(Enum):
    # 1xx - Input / Validation Errors
    MISSING_FIELD = 101             # Required input is missing
    INVALID_TYPE = 102              # Input type mismatch or format error
    INVALID_VALUE = 103             # Input value is invalid or out of range
    DUPLICATE_INPUT = 104           # Redundant/conflicting input data
    INVALID_NAME = 105              # Invalid name format or type

    # 2xx - Authentication / Authorization
    INVALID_TOKEN = 201             # Authentication token is invalid or expired
    UNAUTHORIZED_ACCESS = 202       # User lacks permission for the action
    AUTH_REQUIRED = 203             # Authentication required but not provided

    # 3xx - Resource Errors
    NOT_FOUND = 301                 # Resource does not exist
    ALREADY_EXISTS = 302            # Resource creation conflict
    CONFLICT = 303                  # Resource state conflict (e.g., versioning)

    # 4xx - External / Dependency Failures
    TIMEOUT = 401                   # External service or API timed out
    UPSTREAM_UNAVAILABLE = 402      # Dependent service not reachable
    EXTERNAL_ERROR = 403            # General external service failure

    # 5xx - Internal / System Errors
    INTERNAL_ERROR = 501            # Unhandled exception or unexpected state
    NULL_POINTER = 502              # NoneType access error
    CONFIGURATION_ERROR = 503       # Application misconfigured

    # 6xx - Database / Storage Errors
    DB_CONNECTION_FAILED = 601      # Cannot connect to database
    QUERY_FAILED = 602              # SQL or ORM query failure
    TRANSACTION_ABORTED = 603       # Transaction rollback or failure

    # 7xx - Filesystem / I/O Errors
    FILE_NOT_FOUND = 701            # File or path does not exist
    PERMISSION_DENIED = 702         # File cannot be accessed due to permissions
    WRITE_FAILED = 703              # Write operation failed
    READ_FAILED = 704               # Read operation failed
    FILE_LOAD_FAILED = 705          # Failed to load file

    # 8xx - Business Logic / Domain Errors
    INVALID_STATE = 801             # Operation not allowed in current state
    CONSTRAINT_VIOLATION = 802      # Logical constraint or rule broken
    DUPLICATE_OPERATION = 803       # Duplicate request (e.g., repeated submit)

    # 9xx - Unknown / Other Errors
    UNKNOWN_ERROR = 901             # No known cause, catch-all
    NOT_IMPLEMENTED = 902           # Placeholder for unimplemented feature
    
    def to_dict(self):
        return {"value": self.value, "name": self.name}


class WormcatError(Exception):
    """General-purpose error for Wormcat operations."""

    def __init__(self, message: str, code: dict, origin: str = "", detail: dict = None):
        self.message = message
        self.code_value = code['value']    # int for downstream
        self.code_name = code['name']    # str for logs, JSON, etc.
        self.origin = origin
        if detail is None:
            detail = {}
        self.detail = detail
        
        super().__init__(f"[{self.code_value}] {self.message}")

    def __str__(self):
        origin_info = f" (origin: {self.origin})" if self.origin else ""
        return f"[{self.code_value} {self.code_name}] {self.message} {origin_info}"