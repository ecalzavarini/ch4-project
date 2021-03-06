<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<!-- <xsl:variable name="what">1</xsl:variable> -->


<xsl:template name="treeview">
<xsl:param name="depth"/> 
<xsl:param name="content"/> 
<xsl:param name="radix"/>
<xsl:param name="radix_base"/>

<!-- <xsl:value-of select="number($depth)"/> -->
<xsl:for-each select="$content/core">
<!-- <xsl:value-of select="@name"/> -->
<xsl:variable name="radixc"><xsl:value-of select="$radix"/><xsl:value-of select="@name"/></xsl:variable>
<xsl:choose>
<xsl:when test="number($depth)>=0">
<xsl:if test="number($what)!=3">
#ifdef <xsl:value-of select="$radixc"/>
</xsl:if>
</xsl:when>
</xsl:choose>
<xsl:if test="number($what)=3">
if(<xsl:value-of select="$radixc"/>)
   vmessage(STATUS "<xsl:value-of select="$radixc"/> is ON")
   add_definitions(-D<xsl:value-of select="$radixc"/>)
   <xsl:value-of select="cmake"/>
else(<xsl:value-of select="$radixc"/>)
   vmessage(STATUS "<xsl:value-of select="$radixc"/> is OFF")
endif(<xsl:value-of select="$radixc"/>)
</xsl:if>

<xsl:for-each select="./variable">
<!-- Variable --> <xsl:variable name="radixa"><xsl:value-of select="$radixc"/>_<xsl:value-of select="@name"/></xsl:variable>
<xsl:if test="number($what)=0">

  <xsl:if test="@elements">
    
<!-- Variable --> <xsl:variable name="aaa"><xsl:value-of select="@elements"/></xsl:variable>

  <xsl:if test="@required!='yes'">
    if_defined_get_parameter_<xsl:value-of select="@type"/>_default("<xsl:value-of select="translate($radixa,'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz')"/>",<xsl:text disable-output-escaping="yes">&amp;</xsl:text><xsl:value-of select="translate($radixa,'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz')"/>,<xsl:value-of select="@value"/>,<xsl:value-of select="translate($aaa,'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz')"/>);  </xsl:if>  <xsl:if test="@required='yes'">
    get_parameter_<xsl:value-of select="@type"/>("<xsl:value-of select="translate($radixa,'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz')"/>",<xsl:text disable-output-escaping="yes">&amp;</xsl:text><xsl:value-of select="translate($radixa,'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz')"/>,<xsl:value-of select="translate($aaa,'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz')"/>); </xsl:if>
  </xsl:if>  

  <xsl:if test="not(@elements)">
  <xsl:if test="@required!='yes'">
    if_defined_get_parameter_<xsl:value-of select="@type"/>_default("<xsl:value-of select="translate($radixa,'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz')"/>",<xsl:text disable-output-escaping="yes">&amp;</xsl:text><xsl:value-of select="translate($radixa,'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz')"/>,<xsl:value-of select="@value"/>);  </xsl:if>  <xsl:if test="@required='yes'">
    get_parameter_<xsl:value-of select="@type"/>("<xsl:value-of select="translate($radixa,'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz')"/>",<xsl:text disable-output-escaping="yes">&amp;</xsl:text><xsl:value-of select="translate($radixa,'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz')"/>); </xsl:if>
</xsl:if>



<!-- <xsl:value-of select="$radixa"/> -->
</xsl:if>
<xsl:if test="number($what)=2">
extern <xsl:value-of select="@type"/><xsl:text disable-output-escaping="yes"> </xsl:text><xsl:value-of select="translate($radixa,'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz')"/>;</xsl:if>
<xsl:if test="number($what)=1"><xsl:text>
</xsl:text>
<xsl:value-of select="@type"/><xsl:text disable-output-escaping="yes"> </xsl:text><xsl:value-of select="translate($radixa,'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz')"/>;</xsl:if>
</xsl:for-each>
<xsl:choose>

<xsl:when test="number($depth)=0">
<!-- <xsl:value-of select="@name"/> -->
<xsl:variable name="radixa"><xsl:value-of select="@name"/></xsl:variable>
</xsl:when> 
<xsl:otherwise>
<!--  <xsl:value-of select="$radix"/>_<xsl:value-of select="@name"/> -->
<xsl:variable name="radixa"><xsl:value-of select="$radix"/>_<xsl:value-of select="@name"/></xsl:variable>
</xsl:otherwise>
</xsl:choose>

<xsl:variable name="radixa"><xsl:value-of select="$radix"/><xsl:value-of select="@name"/>_</xsl:variable>
<xsl:variable name="radixb"><xsl:value-of select="$radix"/><xsl:value-of select="@name"/></xsl:variable>
<!-- <xsl:value-of select="$radixa"/> -->
<xsl:call-template name="treeview">
<xsl:with-param name="content" select="."/>
<xsl:with-param name="depth" select="number($depth)+1"/> 
<xsl:with-param name="radix" select="$radixa"/>
<xsl:with-param name="radix_base" select="$radixb"/>
</xsl:call-template>
<xsl:choose>
<xsl:when test="number($depth)>=0">
<xsl:if test="number($what)!=3">
#endif /* <xsl:value-of select="$radixc"/> */
</xsl:if>
</xsl:when>
</xsl:choose>
</xsl:for-each>
</xsl:template>






<xsl:template match="/">

<xsl:if test="number($what)=0">
void read_parameter()
{
  
    read_parameter_file();
</xsl:if>
<xsl:if test="number($what)=2">
#ifndef _COMMON_XML_H_
#define _COMMON_XML_H_
#define vector_double double *
</xsl:if>
<xsl:if test="number($what)=1">
#ifndef _GVARS_XML_H_
#define _GVARS_XML_H_
#define vector_double double *
</xsl:if>
<xsl:call-template name="treeview">
<xsl:with-param name="content" select="/ftmake"/>
<xsl:with-param name="depth" select="0"/>
</xsl:call-template>
<xsl:if test="number($what)=0">
     if (AMIROOT) fclose(fout);
}</xsl:if><xsl:if test="number($what)=2">
#endif /* _COMMON_XML_H_ */
</xsl:if>
<xsl:if test="number($what)=1">
#endif /* _GVARS_XML_H_ */
</xsl:if>
</xsl:template>
</xsl:stylesheet>
