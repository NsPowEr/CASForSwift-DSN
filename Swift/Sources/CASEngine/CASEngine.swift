import Foundation
import CCASEngine

/// Errori emessi dal motore CAS.
public enum CASError: Error, CustomStringConvertible {
    case engineError(message: String)
    case unknownError

    public var description: String {
        switch self {
        case .engineError(let message):
            return "CAS Engine Error: \(message)"
        case .unknownError:
            return "CAS Engine Error: Unknown error"
        }
    }
}

/// Rappresenta un contesto di calcolo CAS.
/// Gestisce l'Arena di memoria (C++) e le definizioni delle variabili.
public class CASContext {
    internal let handle: CASContextRef

    public init() {
        guard let ctx = cas_context_create() else {
            fatalError("Impossibile creare il contesto CAS")
        }
        self.handle = ctx
    }

    deinit {
        cas_context_destroy(handle)
    }

    /// Parsa una stringa in un'espressione simbolica.
    public func parse(_ input: String) throws -> Expression {
        var exprRef: CASExprRef? = nil
        let errorRef = input.withCString { cString in
            cas_parse(self.handle, cString, &exprRef)
        }

        try checkError(errorRef)

        guard let ref = exprRef else {
            throw CASError.unknownError
        }

        return Expression(handle: ref, context: self)
    }

    /// Semplifica un'espressione.
    public func simplify(_ expression: Expression) throws -> Expression {
        var exprRef: CASExprRef? = nil
        let errorRef = cas_simplify(self.handle, expression.handle, &exprRef)

        try checkError(errorRef)

        guard let ref = exprRef else {
            throw CASError.unknownError
        }

        return Expression(handle: ref, context: self)
    }

    /// Valuta numericamente un'espressione (restituisce un'espressione contenente un DecimalLit).
    public func evaluateNumeric(_ expression: Expression) throws -> Expression {
        var exprRef: CASExprRef? = nil
        let errorRef = cas_evaluate_numeric(self.handle, expression.handle, &exprRef)

        try checkError(errorRef)

        guard let ref = exprRef else {
            throw CASError.unknownError
        }

        return Expression(handle: ref, context: self)
    }

    private func checkError(_ errorRef: CASErrorRef?) throws {
        if let err = errorRef {
            let message = String(cString: cas_error_get_message(err))
            cas_error_destroy(err)
            throw CASError.engineError(message: message)
        }
    }
}

/// Rappresenta un'espressione simbolica nel motore CAS.
/// Mantiene un riferimento forte al CASContext per garantire che l'Arena C++ rimanga valida.
public class Expression {
    internal let handle: CASExprRef
    /// Riferimento al contesto per garantire che l'Arena C++ rimanga in vita finché l'espressione è attiva.
    public let context: CASContext

    internal init(handle: CASExprRef, context: CASContext) {
        self.handle = handle
        self.context = context
    }

    deinit {
        cas_expr_destroy(handle)
    }

    /// Restituisce la rappresentazione LaTeX dell'espressione.
    public var latex: String {
        var textPtr: UnsafeMutablePointer<Int8>? = nil
        let errorRef = cas_format_latex(context.handle, handle, &textPtr)
        
        if let err = errorRef {
            let message = String(cString: cas_error_get_message(err))
            cas_error_destroy(err)
            return "Error: \(message)"
        }
        
        if let ptr = textPtr {
            let s = String(cString: ptr)
            cas_string_destroy(ptr)
            return s
        }
        return ""
    }

    /// Restituisce la rappresentazione testuale (parser-friendly) dell'espressione.
    public var text: String {
        var textPtr: UnsafeMutablePointer<Int8>? = nil
        let errorRef = cas_format_text(context.handle, handle, &textPtr)
        
        if let err = errorRef {
            let message = String(cString: cas_error_get_message(err))
            cas_error_destroy(err)
            return "Error: \(message)"
        }
        
        if let ptr = textPtr {
            let s = String(cString: ptr)
            cas_string_destroy(ptr)
            return s
        }
        return ""
    }
}
