/**********************************************************************
Copyright (c) 2013 by Matt Swain <m.swain@me.com>

The MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

***********************************************************************/

#ifndef PYNE_IS_AMALGAMATED
  #include "json.h"
  #include "jsoncustomwriter.h"
#endif

namespace Json {

CustomWriter::CustomWriter( std::string opencurly,
                            std::string closecurly,
                            std::string opensquare,
                            std::string closesquare,
                            std::string colon,
                            std::string comma,
                            std::string indent,
                            int maxWidth)
   : opencurly_( opencurly )
   , closecurly_( closecurly )
   , opensquare_( opensquare )
   , closesquare_( closesquare )
   , colon_( colon )
   , comma_( comma )
   , indent_( indent )
   , maxWidth_( maxWidth )
{
}


std::string
CustomWriter::write( const Value &root )
{
   document_ = "";
   indentString_ = "";
   writeValue( root, document_, false );
   document_ += "\n";
   return document_;
}


void
CustomWriter::writeValue( const Value &value, std::string &doc, bool forceSingleLine )
{
   switch ( value.type() )
   {
   case nullValue:
      doc += "null";
      break;
   case intValue:
      doc += valueToString( value.asLargestInt() );
      break;
   case uintValue:
      doc += valueToString( value.asLargestUInt() );
      break;
   case realValue:
      doc += valueToString( value.asDouble() );
      break;
   case stringValue:
      doc += valueToQuotedString( value.asCString() );
      break;
   case booleanValue:
      doc += valueToString( value.asBool() );
      break;
   case arrayValue:
      {
         bool isMulti = false;
         if (!forceSingleLine)
         {
            std::string valLine = "";
            writeValue( value, valLine, true);
            if (valLine.length() > maxWidth_)
            {
               isMulti = true;
            }
            else
            {
               doc += valLine;
               break;
            }
         }
         doc += opensquare_;
         if (isMulti)
            indent();
         for ( int index =0; index < value.size(); ++index )
         {
            if (isMulti)
            {
               doc += "\n";
               doc += indentString_;
            }
            writeValue( value[index], doc, false );
            if ( index < value.size()-1 )
               doc += comma_;
         }
         if (isMulti)
         {
            unindent();
            doc += "\n";
            doc += indentString_;
         }
         doc += closesquare_;
      }
      break;
   case objectValue:
      {
         bool isMulti = false;
         if (!forceSingleLine)
         {
            std::string valLine = "";
            writeValue( value, valLine, true);
            if (valLine.length() > maxWidth_)
            {
               isMulti = true;
            }
            else
            {
               doc += valLine;
               break;
            }
         }
         Value::Members members( value.getMemberNames() );
         doc += opencurly_;
         if (isMulti)
            indent();
         for ( Value::Members::iterator it = members.begin();
               it != members.end();
               ++it )
         {
            if (isMulti)
            {
               doc += "\n";
               doc += indentString_;

            }
            const std::string &name = *it;
            doc += valueToQuotedString( name.c_str() );
            doc += colon_;
            writeValue( value[name], doc, forceSingleLine );
            if ( !(it + 1 == members.end()) )
               doc += comma_;
         }
         if (isMulti)
         {
            unindent();
            doc += "\n";
            doc += indentString_;
         }
         doc += closecurly_;
      }
      break;
   }
}


void
CustomWriter::indent()
{
   indentString_ += indent_;
}


void
CustomWriter::unindent()
{
   int idSize = int(indent_.size());
   int idsSize = int(indentString_.size());
   if (idsSize >= idSize)
      indentString_.resize (idsSize - idSize);
}

}
