#include "vn/error.h"

#include <string.h>

void strFromVnError(char* out, VnError val)
{
	#if defined(_MSC_VER)
		/* Disable warnings regarding using strcpy_s since this
		 * function's signature does not provide us with information
		 * about the length of 'out'. */
		#pragma warning(push)
		#pragma warning(disable:4996)
	#endif

	switch (val)
	{
	case E_NONE:
		strcpy(out, "NONE");
		break;
	case E_UNKNOWN:
		strcpy(out, "UNKNOWN");
		break;
	case E_BUFFER_TOO_SMALL:
		strcpy(out, "BUFFER_TOO_SMALL");
		break;
	case E_INVALID_VALUE:
		strcpy(out, "INVALID_VALUE");
		break;
	case E_NOT_IMPLEMENTED:
		strcpy(out, "NOT_IMPLEMENTED");
		break;
	case E_NOT_SUPPORTED:
		strcpy(out, "NOT_SUPPORTED");
		break;
	case E_NOT_FOUND:
		strcpy(out, "NOT_FOUND");
		break;
	case E_TIMEOUT:
		strcpy(out, "TIMEOUT");
		break;
	case E_PERMISSION_DENIED:
		strcpy(out, "PERMISSION_DENIED");
		break;
	case E_INVALID_OPERATION:
		strcpy(out, "INVALID_OPERATION");
		break;
	case E_SIGNALED:
		strcpy(out, "SIGNALED");
		break;
	case E_SENSOR_HARD_FAULT:
		strcpy(out, "SENSOR_HARD_FAULT");
		break;
	case E_SENSOR_SERIAL_BUFFER_OVERFLOW:
		strcpy(out, "SENSOR_SERIAL_BUFFER_OVERFLOW");
		break;
	case E_SENSOR_INVALID_CHECKSUM:
		strcpy(out, "SENSOR_INVALID_CHECKSUM");
		break;
	case E_SENSOR_INVALID_COMMAND:
		strcpy(out, "SENSOR_INVALID_COMMAND");
		break;
	case E_SENSOR_NOT_ENOUGH_PARAMETERS:
		strcpy(out, "SENSOR_NOT_ENOUGH_PARAMETERS");
		break;
	case E_SENSOR_TOO_MANY_PARAMETERS:
		strcpy(out, "SENSOR_TOO_MANY_PARAMETERS");
		break;
	case E_SENSOR_INVALID_PARAMETER:
		strcpy(out, "SENSOR_INVALID_PARAMETER");
		break;
	case E_SENSOR_INVALID_REGISTER:
		strcpy(out, "SENSOR_INVALID_REGISTER");
		break;
	case E_SENSOR_UNAUTHORIZED_ACCESS:
		strcpy(out, "SENSOR_UNAUTHORIZED_ACCESS");
		break;
	case E_SENSOR_WATCHDOG_RESET:
		strcpy(out, "SENSOR_WATCHDOG_RESET");
		break;
	case E_SENSOR_OUTPUT_BUFFER_OVERFLOW:
		strcpy(out, "SENSOR_OUTPUT_BUFFER_OVERFLOW");
		break;
	case E_SENSOR_INSUFFICIENT_BAUD_RATE:
		strcpy(out, "SENSOR_INSUFFICIENT_BAUD_RATE");
		break;
	case E_SENSOR_ERROR_BUFFER_OVERFLOW:
		strcpy(out, "SENSOR_ERROR_BUFFER_OVERFLOW");
		break;
	case E_DATA_NOT_ELLIPTICAL:
		strcpy(out, "DATA_NOT_ELLIPTICAL");
		break;
	default:
		strcpy(out, "UNKNOWN_VALUE");
		break;
	}

	#if defined(_MSC_VER)
		#pragma warning(pop)
	#endif
}
